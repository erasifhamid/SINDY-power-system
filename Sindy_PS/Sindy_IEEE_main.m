%% Danish Rafiq and Asif Hamid  
% File Crerated on 19.2.2022
% Recovers system parameters of power-grid models using Sindy
%   Ref:
%   [1] T. Nishikawa and A. E. Motter, Comparative analysis of existing
%       models for power-grid synchronization, New J. Phys. 17, 015012 (2015).
%
%   [2] A. E. Motter, S. A. Myers, M. Anghel, and T. Nishikawa, Spontaneous
%       synchrony in power-grid networks, Nat. Phys. 9, 191-197 (2013).
%
%   [3] Danish Rafiq, Junaid Farooq and M. A. Bazaz, Synergistic use of intrusive and non-intrusive model order
%       reduction techniques for dynamical power grids (2022)
%
%   [4] Sparse identification of non linear dynamic systems SL Brunton (2016)
%% Initialize the IEEE Power System (Load data)
clear;clc;close all
mpc=case4gs;        % 4gs, 57, 118, 145, 300, 1888, 2736sp, test_system_10gen (see data file)
mpc.ref_freq=60;  % reference frequency
global data 
[data,details]=EN_model(mpc); %EN, SM (H,D,A,K,gamma,omega_R)
n_oc=length(data.H); %No of oscillators (FOM Size =2*n_oc)
dt=0.001;tf=5;
tspan=dt:dt:tf;
%  (Full Order model (converted to first order)
f= @(x) power_func(x);
xdotNL= @(t,x) f(x);
%%%%%%%%%% range of initial conditions
x0list = [0]; % need more data to ID
x0list2 = [0]; % need more data to ID
delta = [];
omega = [];
for i = 1:length(x0list)
    x0= [x0list(i)*ones(n_oc,1); x0list2(i)*ones(n_oc,1)]; 
    [tFOM,xFOM]=ode45(xdotNL,tspan,x0);
    delta = [delta; xFOM(:,1:n_oc)];
    omega = [omega; xFOM(:,n_oc+1:end)];
end
D=data.D;
H=data.H;
A=data.A;
gamma=data.gamma;
K=data.K;
omega_R=data.omega_R;
Mi=-D./(2*H);
Di=omega_R./(2*H).*A;
disp('Snapshots recorded')
%% POD (for large data-sets)
% X=xFOM'; Opts.NoSV = 8;
% [xdotr,x0r,Vpod] = POD(X, xdotNL,x0,Opts);
% [tROM,xROM]=ode45(xdotr,tspan,x0r);
% xPOD=Vpod*xROM'; 
%yPOD=C*xPOD./n_oc;
%  subplot(1,2,1)
%  surf(xFOM); shading interp
%  subplot(1,2,2)
%  surf(xPOD'); shading interp
%  figure()
%  plot(tFOM,xFOM,'k',tFOM, xPOD,'r--')
%plot(tFOM,X(1:3,:),tFOM, xPOD(1:3,:),'k--')
%% Sindy 
N=n_oc;
x=omega;
%create dx using finite-difference scheme
%dx=zeros(size(x));
% for i=1:size(x,2)
%     dx(1,i)=(x(2,i)-x(1,i))/dt;
%     dx(end,i)=(x(end,i)-x(end-1,i))/dt;
%         for j=2:length(dxx)-1
%             dx(j,i)=(x(j+1,i)-x(j-1,i))/(2*dt);
%         end
% end

% another option is to creare using the known nonlinear function
xFOM=[delta omega];
for i=1:length(xFOM)
    dx(i,:)=xdotNL(0,xFOM(i,:)');
end
dx=dx(:,n_oc+1:end);
lambda=0.1; % 0.1
polyorder=1;
Theta=poolDataie(x,polyorder,delta,gamma,K);
Xi = sparsifyDynamics(Theta,dx,lambda,n_oc); 
%% Reconstruction
xall2 = [];
dxall2 = [];
for i = 1:length(x0list)
    x0= x0list(i)*ones(n_oc,1);
    [tB,xB]=ode45(@(t,x)SparseGalerkinie(t,x,tspan,Xi,polyorder,delta,gamma,K),tspan,x0(1:n_oc));% approximate
    xall2 = [xall2; xB];
end
xs=xall2;
%
figure(1);
%f.Position = [200 200 700 600];
%subplot(2,1,1)
plot(dt*(1:length(x)),x,'k', dt*(1:length(xs)),xs,'r--','linewidth',1.5)
xlabel('$t$','Interpreter','latex')
ylabel('$\omega$','Interpreter','latex')
grid on
set(gca,'ticklabelinterpreter','latex','Fontsize',15)
%figure(2)
%subplot(2,1,2)
% plot(tspan,xFOM(:,n_oc+1:end),'g')
% hold on
% plot(tspan,xB(:,n_oc+1:end),'k--','linewidth',1.5)
% xlabel('$t$','Interpreter','latex')
% ylabel('$\omega$','Interpreter','latex')
% grid on
% set(gca,'ticklabelinterpreter','latex','Fontsize',15)

%%
norm_full=norm(x-xs)
D=data.D;
H=data.H;
A=data.A;
gamma=data.gamma;
K=data.K;
omega_R=data.omega_R;
const_coeff = (omega_R./(2*H).*A)'
x_coeff = - diag(D./(2*H))
sin_coeff= - diag(omega_R./(2*H))

