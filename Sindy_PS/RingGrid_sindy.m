%  SinDy test on Ring-Grid Model
% Danish Rafiq
% File created on 01.03.2022
% The script solves Ring Grid Power System Model using SinDy algorithm
% Eqn: delta_dd +0.25*delta_d = 0.95 - sin(delta) -100(sin(delta_i-delta_i+1) + sin(delta_i - delta_i-1))
clear all;clc;close all
% System Parameters
N=1;         % Full Order of the system
u_train=@(t) 0.52; %0.12-0.52
u_test = u_train;
B=ones(2*N,1);
C=zeros(1,2*N);
C(1)=1;
IC=[1.15*ones(1,N) zeros(1,N)]'; %initial conditions 
dt=0.01;
tsim=0:dt:30;    % Simulation time
f=@(x) Swing_equation(x);
xdotNL= @(t,x)f(x)+u_train(t);
disp('Taking Snapshots of the system:...')
[tFOM,xFOM]=ode45(xdotNL,tsim,IC);
%% compute Derivative 
%eps = .05;      % noise strength
x=xFOM;
for i=1:length(x)
    dx(i,:) = xdotNL(0,x(i,:)) ;
end
%% pool Data  (i.e., build library of nonlinear time series)
% polyorder = 3;
% usesine = 1;
Theta = poolDataPS(x,N);
% m = size(Theta,2);
%% compute Sparse regression: sequential least squares
lambda = 0.0001;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,N)
%poolDataLISTNew({'x'},Xi,N);
%poolDataLIST({'x','y'},Xi,2,1,0);
%%
% [tB,xB]=ode45(@(t,x)sparseGalerkinNew(x,Xi,polyorder,usesine),tsim,IC); 
[tB,xB]=ode45(@(t,x)sparseGalerkinPS(x,Xi),tsim,IC);% approximate

%% Plotting
FS=15;INT='latex';
f=figure(1);
f.Position = [200 200 700 600];
subplot(2,1,1)
plot(tsim,xFOM(:,1),'r','linewidth',2)
hold on
plot(tsim,xFOM(:,2),'c','linewidth',2)
plot(tsim,xB(:,1),'k--','linewidth',2)
plot(tsim,xB(:,2),'k--','linewidth',2)
ylim([-2 2])
grid on
l=legend('True $\delta$','True $\omega$','Identified');
ylabel('States','interpreter',INT); xlabel('$t$','interpreter',INT);
set(gca,'ticklabelinterpreter',INT,'Fontsize',FS)
set(l,'interpreter','latex','location','northeast')


subplot(2,1,2)
plot(xFOM(:,1),xFOM(:,2),'g','LineWidth',2)
hold on
plot(xB(:,1),xB(:,2),'k--','Linewidth',2)
xlabel('$\delta$','interpreter','latex');
ylabel('$\omega$','Interpreter','latex'); xlim([-0.6 1.2])
legend('True', 'Identified')
grid on
set(gca,'ticklabelinterpreter',INT,'Fontsize',FS)
%% plotting
% figure
% subplot(1,2,1)
% dtA = [0; diff(tA)];
% color_line3(xA(:,1),xA(:,2),xA(:,3),dtA,'LineWidth',1.5);
% view(27,16)
% grid on
% xlabel('x','FontSize',13)
% ylabel('y','FontSize',13)
% zlabel('z','FontSize',13)
% set(gca,'FontSize',13)
% subplot(1,2,2)
% dtB = [0; diff(tB)];
% color_line3(xB(:,1),xB(:,2),xB(:,3),dtB,'LineWidth',1.5);