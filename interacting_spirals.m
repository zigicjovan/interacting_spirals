% Clean workspace
clear all; close all; clc

%% General spiral wave data using spectral methods

% System parameters
d1 = 0.1; 
d2 = 0.1; 
beta = 1.0;

tic
% Spectral method variables
Tmax = 10;                                                                 % maximum time window to save in one file
T = 100;                                                                   % time window to simulate
dt = 0.05;                                                                 % time-step size    
Ntime = T/dt + 1;                                                          % total number of time states
Ntime_save_max = Tmax/dt ;                                              % maximum length of saved data file
timewindow = linspace(0,T,Ntime);
Tstops = T;
if T > Tmax                                                                % create intervals to compute if larger than Tmax
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
else                                                                       % small windows just need one file
    intervals = 1;
    tspan = 0:dt:T;
end
Lx = 40;                                                                    % size of X-dim
Ly = 40;                                                                    % size of Y-dim 
n = 256;                                                                    % spatial grid resolution (number of Fourier modes in each dimension)
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

Y = Y + 10;
X = X + 10;
u0 = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));         % initial condition in u
v0 = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+1i*Y)-(sqrt(X.^2+Y.^2)));         % initial condition in v
u0 = [fliplr(u0) u0];                                                       % counter-rotation in u    
v0 = [fliplr(v0) v0];                                                       % counter-rotation in u
%u0 = [u0 u0];                                                               % co-rotation in u
%v0 = [v0 v0];                                                               % co-rotation in v
uv0 = [reshape(fft2(u0),1,N) reshape(fft2(v0),1,N)].';                      % initial condition in Fourier space
try
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(100) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
    uv0 = readmatrix(four_file)';
catch
end 

%
% Numerical timestepping in Fourier domain
mkdir([pwd  '/data/forward' ]);
uvI = uv0;
for sim = 1:intervals-1                                                     % do for all but last interval
    [t, uvsol] = ode45(@(t,uvt) rdpde(t,uvt,K22,d1,d2,beta,n,N),[ tspan(sim,:) tspan(sim,end) + dt ],uvI);
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(tspan(sim,end)) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
    uvI = conj(uvsol(end,:)');
    writematrix(uvsol(1:end-1,:), four_file);
end
tspanLast = rmmissing(tspan(end,:));                                       % deal with last interval which may be different length
if length(tspanLast) >  1
    [t, uvsol] = ode45(@(t,uvt) rdpde(t,uvt,K22,d1,d2,beta,n,N),tspanLast,uvI);
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(T) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
    writematrix(uvsol, four_file);
else
    uvsol = uvsol(end,:);
    four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(T) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
    writematrix(uvsol, four_file);
end
%}

%{
fid = fopen('filename.bin','w');
fwrite(fid,var,'double');
fclose(fid);
%}

toc
%{
%% Reshape and view spiral wave
figure;
set(gcf,'Position',[100 100 900 750])
axis tight manual % this ensures that getframe() returns a consistent size
mkdir([pwd  '/media/movies' ]);

filename = [pwd '/media/movies/phys_n_' num2str(n) ...
            '_T_' num2str(T) '_dt_' num2str(d1) '_L1_' num2str(Lx,'%.3f') '_L2_' num2str(Ly,'%.3f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.gif'];

title1 = 'Forward-time Ginzburg-Landau solution';
Ntime = size(uvsol,1);
timewindow = linspace(0,T,Ntime);
set(gcf,'color','white')
set(gca,'color','white')

frames = ceil(Ntime/T);

for j = 1 : frames-1 : length(tspan)
    currentT = timewindow(1,j);
    title2 = ['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(currentT,'%.1f') ', N = ' num2str(n) '$'];

    ut = reshape((uvsol(j,1:N).'),n,n);
    vt = reshape((uvsol(j,(N+1):(2*N)).'),n,n);
    u(:,:,j) = real(ifft2(ut));
    v(:,:,j) = real(ifft2(vt));

    pcolor(u(:,:,j)); 
    %pcolor(x,y,real(ifft2(reshape((uvsol(j,1:N).'),n,n))))
    shading interp; colormap(hot); colorbar; drawnow; 
    title({title1, title2},'Interpreter','latex', 'FontSize',14);
    if j == 1
        gif(filename)
    else
        gif
    end
end
%}

%% post-processing

% compute L2 norm and Fourier mode evolution
normL2 = NaN(Ntime,1);
v_mean = zeros(round(sqrt((n/2)^2+(n/2)^2)) + 1,Ntime);
v_meancount = v_mean;
stopcounter = 1;
currentTstop = Tstops(stopcounter,1);
four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
uvsol = readmatrix(four_file);
for i = 1:Ntime
    
    currentT = timewindow(1,i);
    if currentT > currentTstop && currentTstop ~= Tstops(end,1)
        stopcounter = stopcounter + 1;
        currentTstop = Tstops(stopcounter,1);
        four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
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

mkdir([pwd  '/data' ]);
mkdir([pwd  '/data/normL2' ]);
mkdir([pwd  '/data/spectrum' ]);
normL2data_file = [pwd '/data/normL2/normL2_n_' num2str(n) ...
            '_T_' num2str(T) '_dt_' num2str(d1) '_L1_' num2str(Lx,'%.3f') '_L2_' num2str(Ly,'%.3f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
spectrum_file = [pwd '/data/spectrum/spectrum_n_' num2str(n) ...
            '_T_' num2str(T) '_dt_' num2str(d1) '_L1_' num2str(Lx,'%.3f') '_L2_' num2str(Ly,'%.3f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];

writematrix(normL2, normL2data_file);
writematrix(v_mean, spectrum_file);

%
%% make gif
figure;
set(gcf,'Position',[100 100 1800 1500])
axis tight manual % this ensures that getframe() returns a consistent size
mkdir([pwd  '/media/movies' ]);

filename = [pwd '/media/movies/diagnostics_n_' num2str(n) ...
            '_T_' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_frames_' num2str(length(timewindow)) '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.gif'];

title1 = 'Forward-time Ginzburg-Landau solution';
set(gcf,'color','white')
set(gca,'color','white')

frames = ceil(Ntime/T);
stopcounter = 1;
framecounter = 1;
currentTstop = Tstops(stopcounter,1);
four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
uvsol = readmatrix(four_file);
for i = 1 : frames-1 : Ntime

    currentT = timewindow(1,i);
    if currentT > currentTstop && currentTstop ~= Tstops(end,1)
        stopcounter = stopcounter + 1;
        currentTstop = Tstops(stopcounter,1);
        four_file = [pwd '/data/forward/four_n_' num2str(n) '_T_' num2str(currentTstop) '_Lx_' num2str(Lx,'%.0f') '_Ly_' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') '.dat'];
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
    pcolor(u(:,:,i)); 
    %pcolor(x,y,real(ifft2(reshape((uvsol(j,1:N).'),n,n))))
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
    %legend('L^{2} norm','Location','southeast')

    subplot(2,3,6);
    semilogy(v_mean(:,i),".")
    xlabel('$k \approx \sqrt{k_1^2+k^2_2}$','Interpreter','latex','FontSize',12); 
    ylabel('$\frac{1}{j}\sum_{j} |{\widehat{u}_k}|$','Interpreter','latex','FontSize',12);
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

mkdir([pwd  '/media/energy' ]);
% Wavenumber evolution plot
timewindow = linspace(0,T,Ntime);
h = figure;
semilogy(timewindow,v_mean(1,:),'LineWidth',0.1,'Marker','.')
hold on;
for i = 2:size(v_mean,1)
    semilogy(timewindow,v_mean(i,:),'LineWidth',0.1,'Marker','.')
end
set(gcf,'Position',[100 100 900 750])
xlabel('Time $t$','Interpreter','latex'); 
xlim([0 T])
ylim([1e-15 max(v_mean(1,:))+1e5 ])
ylabel('$\frac{1}{j}\sum_{j} |{\widehat{u}_k}|$','Interpreter','latex');
fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title('Evolution of Fourier spectrum','Interpreter','latex')
subtitle(['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(T,'%.1f') ', N = ' num2str(n) '$'],'Interpreter','latex','FontSize',14)
filename = [pwd '/media/energy/wavenumberevol_n' num2str(n) '_T' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ''];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

% L2 norm plot
h = figure;
semilogy(timewindow,normL2,'LineWidth',0.5,'Marker','.')
set(gcf,'Position',[100 100 900 750])
xlabel('Time $t$','Interpreter','latex'); 
xlim([0 T])
ylabel('$||{u(t;\varphi)}||_{L^2}$','Interpreter','latex');
fontsize(12,"points")
set(gca,'fontsize', 16) 
set(gcf,'color','white')
set(gca,'color','white')    
title('Evolution of $L^2$ norm','Interpreter','latex')
subtitle(['$L_x = ' num2str(Lx,'%.0f') ', L_y = ' num2str(Ly,'%.0f') ', T = ' num2str(T,'%.1f') ', N = ' num2str(n) '$'],'Interpreter','latex','FontSize',14)
filename = [pwd '/media/energy/normL2_n' num2str(n) '_T' num2str(T) '_Lx' num2str(Lx,'%.0f') '_Ly' num2str(Ly,'%.0f') '_d1_' num2str(d1,'%.1f') '_d2_' num2str(d2,'%.1f') ''];
saveas(h,[filename '.fig'])
exportgraphics(h,[filename '.pdf'])

%% Reaction-diffusion RHS
function rhs = rdpde(t,uvt,K22,d1,d2,beta,n,N)

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