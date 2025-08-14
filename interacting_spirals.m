% Clean workspace
clear all; close all; 

if not(isfolder([pwd  '/data' ]))
    mkdir([pwd  '/data' ])
    mkdir([pwd  '/data/forward' ])
    mkdir([pwd  '/data/normL2' ])
    mkdir([pwd  '/data/spectrum' ]);
    mkdir([pwd  '/media' ]);
    mkdir([pwd  '/media/movies' ]);
    mkdir([pwd  '/media/energy' ]);
end

%% General spiral wave data using spectral methods

% System parameters
d1 = 0.1; 
d2 = 0.1; 
beta = 1.0;
continuefrom = 0;                                                           % choose if IC is continued from file
transX = 0.00; transY = 0.00;                                               % translate spiral foci
rotate = 'contra';                                                          % choose spiral rotation: 'co' or 'contra' 
flip = 'no';                                                                % choose spiral flip: 'yes' or 'no'
pert = 1; pert_amp = 1/5; x_amp = 1; y_amp = 1;                             % adjust perturbation magnitude
shiftxr = 0.15; shiftyr = 0.23; shiftxl = 0.00; shiftyl = 0.00;             % translate perturbation foci

tic
% Spectral method variables
Tmax = 10;                                                                  % maximum time window to save in one file
T = 200;                                                                    % time window to simulate
dt = 0.05;                                                                  % time-step size    
Ntime = T/dt + 1;                                                           % total number of time states
Ntime_save_max = Tmax/dt ;                                                  % maximum length of saved data file
timewindow = linspace(0,T,Ntime);

% create intervals to compute solution
current0 = 0;
Tremaining = T;
intervals = ceil(T/Tmax + 1e-10);
Tstops = NaN(intervals,1);
tspan = NaN(intervals,Tmax/dt);
for  i = 1:intervals-1
    Tremaining = Tremaining - Tmax;
    tspan(i,:) = current0:dt:(T-Tremaining-dt);
    current0 = current0+Tmax;
    Tstops(i,1) = current0 - dt;
end
tspan(intervals,1:(T-current0)/dt+1) = current0:dt:T;
Tstops(end,1) = T;

% problem and grid parameters
Lx = 30;                                                                    % size of X-dim
Ly = 30;                                                                    % size of Y-dim 
n = 2^8;                                                                    % spatial grid resolution (number of Fourier modes in each dimension)
N = n*n;                                                                    % total number of grid points
x2 = linspace(-Lx/2,Lx/2,n+1); 
x = x2(1:2:n); 
y2 = linspace(-Ly/2,Ly/2,n+1); 
y = y2(1:n); 
kx = (2*pi/Lx)*[0:(n/2-1) -n/2:-1]; 
ky = (2*pi/Ly)*[0:(n/2-1) -n/2:-1]; 

% Initial Conditions
m = 1;                                                                      % number of spirals
[X,Y] = meshgrid(x,y);
[KX,KY] = meshgrid(kx,ky);
K2 = KX.^2 + KY.^2; 
K22=reshape(K2,N,1);
u = zeros(length(kx),length(ky),length(tspan));
v = zeros(length(kx),length(ky),length(tspan));

if continuefrom > 0
    %four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(continuefrom) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(continuefrom) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
    uv0 = readmatrix(four_file)';
else
    Y = Y + transX*Ly;                                                      % add translation in X
    X = X + transY*Lx;                                                      % add translation in Y
    u0r = tanh(sqrt(X.^2+Y.^2)) .* cos(m*angle(X + 1i*Y) - sqrt(X.^2+Y.^2)) ...
        + pert*exp(-pert_amp*( (x_amp*(X-shiftxr*Lx)).^2 + (y_amp*(Y-shiftyr*Ly)).^2)); % initial condition in u (right side)
    u0l = tanh(sqrt(X.^2+Y.^2)) .* cos(m*angle(X + 1i*Y) - sqrt(X.^2+Y.^2)) ...
        + pert*exp(-pert_amp*( (x_amp*(X-shiftxl*Lx)).^2 + (y_amp*(Y-shiftyl*Ly)).^2)); % initial condition in u (left side)
    v0r = tanh(sqrt(X.^2+Y.^2)) .* sin(m*angle(X + 1i*Y) - sqrt(X.^2+Y.^2)) ...
        + pert*exp(-pert_amp*( (x_amp*(X-shiftxr*Lx)).^2 + (y_amp*(Y-shiftyr*Ly)).^2)); % initial condition in v (right side)
    v0l = tanh(sqrt(X.^2+Y.^2)) .* sin(m*angle(X + 1i*Y) - sqrt(X.^2+Y.^2)) ...
        + pert*exp(-pert_amp*( (x_amp*(X-shiftxl*Lx)).^2 + (y_amp*(Y-shiftyl*Ly)).^2)); % initial condition in v (left side)
    switch rotate
        case 'co'
            u0 = [u0l u0r];                                                 % co-rotation in u    
            v0 = [v0l v0r];                                                 % co-rotation in v
        case 'contra'
            u0 = [fliplr(u0l) u0r];                                         % contra-rotation in u    
            v0 = [fliplr(v0l) v0r];                                         % contra-rotation in v
    end
    switch flip
        case 'yes'
            switch 'rotate'
                case 'co'
                    u0 = [flipud(u0l) u0r];                                 % co-rotation and flip in u    
                    v0 = [flipud(v0l) v0r];                                 % co-rotation and flip in v
                case 'contra'
                    u0 = [rot90(u0l,2) u0r];                                % contra-rotation and flip in u    
                    v0 = [rot90(v0l,2) v0r];                                % contra-rotation and flip in v
            end
    end
    uv0 = [reshape(fft2(u0),1,N) reshape(fft2(v0),1,N)].';                  % initial condition in Fourier space
end 

%
% Numerical timestepping in Fourier domain
uvI = uv0;
toclist = NaN(intervals,1);
for sim = 1:intervals-1                                                     % do for all but last interval
    [t, uvsol] = ode45(@(t,uvt) rdpde(t,uvt,K22,d1,d2,beta,n,N),[ tspan(sim,:) tspan(sim,end) + dt ],uvI);
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(tspan(sim,end)) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
    uvI = conj(uvsol(end,:)');
    writematrix(uvsol(1:end-1,:), four_file);
    toclist(sim,1) = toc;
    toclist(sim,1) = toclist(sim,1) - sum(toclist(1:sim-1,1));
    remainingtime = toclist(2,1)*(intervals-sim-1);
    if isnan(toclist(2,1))
        remainingtime = 0.9*toclist(1,1)*(intervals-sim-1);
    end
    disp(['Solved up to T = ' num2str(Tstops(sim) + dt) ' after ' num2str(sum(toclist(1:sim,1)),'%.2f') ' s, estimated remaining time is ' num2str(remainingtime,'%.2f') ' s...']);
end
tspanLast = rmmissing(tspan(end,:));                                       % deal with last interval which may be different length
toclist(intervals) = toc;
if length(tspanLast) >  1
    [t, uvsol] = ode45(@(t,uvt) rdpde(t,uvt,K22,d1,d2,beta,n,N),tspanLast,uvI);
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(T) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
    writematrix(uvsol, four_file);
    disp(['Completed solution for T = ' num2str(Tstops(intervals)) ' after ' num2str(toclist(intervals),'%.2f') ' s.']);
else
    uvsol = uvsol(end,:);
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(T) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
    writematrix(uvsol, four_file);
    disp(['Completed solving to T = ' num2str(Tstops(intervals)) ' after ' num2str(toclist(intervals),'%.2f') ' s.']);
end
%}

%{
fid = fopen('filename.bin','w');
fwrite(fid,var,'double');
fclose(fid);
%}

%% post-processing
tic
disp(['Post-processing, estimated remaining time is ' num2str(1.05*toclist(intervals),'%.2f') ' s.']);

% initial state plot
four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(Tstops(1)) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
uvsol = readmatrix(four_file);
h = figure;
set(gcf,'Position',[100 100 900 750])
set(gcf,'color','white')
set(gca,'color','white')  
ut = reshape((uvsol(1,1:N).'),n,n);
pcolor(linspace(0,Lx,n),linspace(0,Ly,n),real(ifft2(ut))); 
xlabel('$x$','Interpreter','latex','FontSize',18);
ylabel('$y$','Interpreter','latex','FontSize',18);
shading interp; colormap(hot); colorbar; drawnow; 
title1 = 'Initial state of Ginzburg-Landau solution';
title2 = ['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(0,'%.0f') ', N = ' num2str(n) '$'];
title({title1, title2},'Interpreter','latex','FontSize',18);
filename = [pwd '/media/energy/initial_n' num2str(n) '_T' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') ''];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

% terminal state plot
four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(Tstops(end)) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
uvsol = readmatrix(four_file);
h = figure;
set(gcf,'Position',[100 100 900 750])
set(gcf,'color','white')
set(gca,'color','white')  
ut = reshape((uvsol(end,1:N).'),n,n);
pcolor(linspace(0,Lx,n),linspace(0,Ly,n),real(ifft2(ut))); 
xlabel('$x$','Interpreter','latex','FontSize',18);
ylabel('$y$','Interpreter','latex','FontSize',18);
shading interp; colormap(hot); colorbar; drawnow; 
title1 = 'Terminal state of Ginzburg-Landau solution';
title2 = ['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(T,'%.0f') ', N = ' num2str(n) '$'];
title({title1, title2},'Interpreter','latex','FontSize',18);
filename = [pwd '/media/energy/terminal_n' num2str(n) '_T' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') ''];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

% compute L2 norm and Fourier mode evolution
normL2 = NaN(Ntime,1);
v_mean = zeros(round(sqrt((n/2)^2+(n/2)^2)) + 1,Ntime);
v_meancount = v_mean;
stopcounter = 1;
currentTstop = Tstops(stopcounter,1);
four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
uvsol = readmatrix(four_file);
for i = 1:Ntime
    
    currentT = timewindow(1,i);
    if currentT > currentTstop && currentTstop ~= Tstops(end,1)
        stopcounter = stopcounter + 1;
        currentTstop = Tstops(stopcounter,1);
        four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
        uvsol = readmatrix(four_file);
    end

    if mod(i,Ntime_save_max) ~= 0 
        imod = mod(i,Ntime_save_max);
    else
        imod = Ntime_save_max;
    end

    ut = reshape((uvsol(imod,1:N).'),n,n);
    u(:,:,i) = real(ifft2(ut));
    u_i = u(:,:,i);
    normL2(i,1) = sqrt(sum( u_i(:) .* conj(u_i(:)) )*(Lx*Ly)/(N*2));
    v = fftshift(real(abs(fft2(u_i))));
    for j = 1:n
        for k = 1:n
            index = round(sqrt((j-(n/2+1))^2+(k-(n/2+1))^2)) + 1;
            v_mean(index,i) = v_mean(index,i) + v(j,k);
            v_meancount(index,i) = v_meancount(index,i) + 1;
        end
    end
    for m = 1:size(v_meancount,1)
        v_mean(m,i) = v_mean(m,i)/v_meancount(m,i);
    end
        
    if i == 1
        u_IC = u_i;
    elseif i == Ntime
        u_TC = u_i;
    end

end
v_mean = v_mean(2:end,:);

normL2data_file = [pwd '/data/normL2/normL2_n_' num2str(n) ...
            '_T_' num2str(T) '_dt_' num2str(d1) '_L1_' num2str(Lx,'%.3f') '_L2_' num2str(Ly,'%.3f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
spectrum_file = [pwd '/data/spectrum/spectrum_n_' num2str(n) ...
            '_T_' num2str(T) '_dt_' num2str(d1) '_L1_' num2str(Lx,'%.3f') '_L2_' num2str(Ly,'%.3f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];

writematrix(normL2, normL2data_file);
writematrix(v_mean, spectrum_file);

% Wavenumber evolution plot
timewindow = linspace(0,T,Ntime);
h = figure;
semilogy(timewindow,v_mean(1,:),'LineWidth',0.1,'Marker','.')
hold on;
for i = 2:size(v_mean,1)
    semilogy(timewindow,v_mean(i,:),'LineWidth',0.1,'Marker','.')
end
set(gcf,'Position',[100 100 900 750])
xlabel('Time $t$','Interpreter','latex','FontSize',18); 
xlim([0 T])
ylim([1e-15 max(v_mean(1,:))+1e5 ])
ylabel('$\frac{1}{j}\sum_{j} |{\widehat{u}_k}|$','Interpreter','latex','FontSize',18);
set(gcf,'color','white')
set(gca,'color','white')    
title('Evolution of Fourier spectrum','Interpreter','latex','FontSize',18)
subtitle(['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(T,'%.1f') ', N = ' num2str(n) '$'],'Interpreter','latex','FontSize',14)
filename = [pwd '/media/energy/wavenumberevol_n' num2str(n) '_T' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') ''];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

% L2 norm plot
h = figure;
semilogy(timewindow,normL2,'LineWidth',0.5,'Marker','.')
set(gcf,'Position',[100 100 900 750])
xlabel('Time $t$','Interpreter','latex','FontSize',18); 
xlim([0 T])
ylabel('$||{u(t;\varphi)}||_{L^2}$','Interpreter','latex','FontSize',18);
set(gcf,'color','white')
set(gca,'color','white')    
title('Evolution of $L^2$ norm','Interpreter','latex','FontSize',18)
subtitle(['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(T,'%.1f') ', N = ' num2str(n) '$'],'Interpreter','latex','FontSize',14)
filename = [pwd '/media/energy/normL2_n' num2str(n) '_T' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') ''];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%
% make gif
figure;
set(gcf,'Position',[100 100 1800 1500])
axis tight manual % this ensures that getframe() returns a consistent size

filename = [pwd '/media/movies/diagnostics_n_' num2str(n) ...
            '_T_' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.gif'];

title1 = 'Forward-time Ginzburg-Landau solution';
set(gcf,'color','white')
set(gca,'color','white')

frames = ceil(Ntime/T);
stopcounter = 1;
framecounter = 1;
currentTstop = Tstops(stopcounter,1);
four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
uvsol = readmatrix(four_file);
for i = 1 : frames-1 : Ntime

    currentT = timewindow(1,i);
    if currentT > currentTstop && currentTstop ~= Tstops(end,1)
        stopcounter = stopcounter + 1;
        currentTstop = Tstops(stopcounter,1);
        four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ...
        '_continuefrom_' num2str(continuefrom,'%.0f') '_transX_' num2str(transX,'%.2f') '_transY_' num2str(transY,'%.2f') '_rot_' rotate '_flip_' flip '_pert_' num2str(pert,'%.0f') '_pertamp_' num2str(pert_amp,'%.1f') '_xamp_' num2str(x_amp,'%.1f') '_yamp_' num2str(y_amp,'%.1f') ...
        '_sxr_' num2str(shiftxr,'%.2f') '_syr_' num2str(shiftyr,'%.2f') '_sxl_' num2str(shiftxl,'%.2f') '_syl_' num2str(shiftyl,'%.2f') '.dat'];
        uvsol = readmatrix(four_file);
    end

    if mod(i,Ntime_save_max) ~= 0 
        imod = mod(i,Ntime_save_max);
    else
        imod = Ntime_save_max;
    end

    title2 = ['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(currentT,'%.1f') ', N = ' num2str(n) '$'];

    subplot(2,3,[1,2,4,5]);
    ut = reshape((uvsol(imod,1:N).'),n,n);
    vt = reshape((uvsol(imod,(N+1):(2*N)).'),n,n);
    u(:,:,i) = real(ifft2(ut));
    v(:,:,i) = real(ifft2(vt));
    pcolor(linspace(0,Lx,n),linspace(0,Ly,n),u(:,:,i)); 
    xlabel('$x$','Interpreter','latex','FontSize',18);
    ylabel('$y$','Interpreter','latex','FontSize',18);
    shading interp; colormap(hot); colorbar; drawnow; 
    title({title1, title2},'Interpreter','latex','FontSize',18);

    subplot(2,3,3);
    semilogy(timewindow,normL2,'b')
    hold on
    xline(currentT,'-');
    hold off
    xlabel('Time $t$','Interpreter','latex','FontSize',18);
    ylabel('$|| u(t;\varphi) ||$','Interpreter','latex','FontSize',18);
    xlim([0 T])
    ylim([min(normL2) max(normL2)+1e-1])
    title("Evolution of $L^2$ norm",'Interpreter','latex','FontSize',18)

    subplot(2,3,6);
    semilogy(v_mean(:,i),".")
    xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex','FontSize',18); 
    ylabel('$\frac{1}{j}\sum_{j} |{\widehat{u}_k}|$','Interpreter','latex','FontSize',18);
    title("Energy spectrum",'Interpreter','latex','FontSize',18)
    xlim([1 size(v_mean,1)])
    ylim([1e-15 max(max(v_mean))])
    
    if i == 1
        gif(filename)
    else
        gif
    end
end
%}
disp(['Post-processing completed after ' num2str(toc,'%.2f') ' s.']);

%% Reaction-diffusion RHS
function rhs = rdpde(~,uvt,K22,d1,d2,beta,n,N)

    % Calculate u and v terms
    ut = reshape((uvt(1:N)),n,n);
    vt = reshape((uvt((N+1):(2*N))),n,n);
    u = real(ifft2(ut)); 
    v = real(ifft2(vt));

    % Reaction Terms
    u3 = u.^3; 
    v3 = v.^3; 
    u2v = (u.^2).*v; 
    uv2 = u.*(v.^2);
    utrhs = reshape((fft2(u-u3-uv2+beta*u2v+beta*v3)),N,1);
    vtrhs = reshape((fft2(v-u2v-v3-beta*u3-beta*uv2)),N,1);

    rhs = [-d1*K22.*uvt(1:N)+utrhs
         -d2*K22.*uvt(N+1:end)+vtrhs];
 
end