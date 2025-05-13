% This file is a part of manuscript "Structure Preserving Nonlinear Reduced Order Modeling Technique for Power Systems"
% Authors: Danish Rafiq and Mohammad Abid Bazaz
% Contact: danish_pha2007@nitsri.ac.in

% Solves the Ring Grid Power System Model for three scenarios (as discussed in the paper) in second-order structure
% with n similar generators via DMD & POD using implicit Euler Method
% M*delta_ddot + D* delta_dot = pm-(b*sin(delta_i)+bint[sin(delta_i - delta_i+1)+sin(delta_i - delta_i-1)])

% Selected References:
% [1] Reduced order modelling for transient simulation of power systems
%     using TPWL (Haris Malik, Borzacchiello,Francisco and Diez)
% [2] Nonlinear Moment Matching for the Simulation-Free Reduction of
%     Structural Systems (Maria et al.)

%%
clear all; close all; clc
n=100;      % Dimension of Full Order Model
rdefl=15;    % Dimension of Reduced Model
u=@(t)0.95;  % Pm in pu 
tf=50;       
del=1e-3;       
tSim=0:del:tf-del;  % Simulation time
k=round(tf/del);
%FOM
M=eye(n);
D=0.25*eye(n);
B=ones(n,1);
C=ones(1,n);
f=@(x)-sin(x) - 100*(sin(x-circshift(x,1)) + sin(x-circshift(x,-1)));
Jf=@(x) -cos(x) - 100*(cos(x-circshift(x,1)) + cos(x-circshift(x,-1)));
% Case1: All nodes starting from under-perturbation
x0=0.8*ones(n,1); % initial conditions
steps=k;
disp('Simulating Case-I: All nodes starting from under perturbation..')
[xFOM1,FOM_time_case1]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM1=C*xFOM1./n;  
%%
%% POD
[Vpod,SS,VV]=svd(xFOM1,'econ');
Phi=Vpod(:,1:4);
Mr=Phi'*M*Phi;
Dr=Phi'*D*Phi;
Br=Phi'*B;
Cr=C*Phi;
xr=Phi'*x0;
x=@(x)Phi*x;
%fr=@(x)-sin(x) - 100*(sin(x-circshift(x,1)) + sin(x-circshift(x,-1)))
[xROM_pod,ROM_pod_time1]=implicitEulerSO(Mr,Dr,Br,xr,f,del,u,steps);
xROM=xROM_pod';
%yROM_pod1=C*xROM_pod./n;

%% DMD
% X = xFOM1(:,1:end-1);
% X2 = xFOM1(:,2:end);
% [U,S,V] = svd(X,'econ');
% %  Compute DMD (Phi are eigenvectors)
% r = 150;  % truncate at 15 modes
% U = U(:,1:r);
% S = S(1:r,1:r);
% V = V(:,1:r);
% Atilde = U'*X2*V*inv(S);
% [W,eigs] = eig(Atilde);
% Phi = X2*V*inv(S)*W;
% [Phi,~] = GramSchmidt(Phi);
% %% ROM using DMD
% Mr=Phi'*M*Phi;
% Dr=Phi'*D*Phi;
% Br=Phi'*B;
% Cr=C*Phi;
% xn=ones(n,1);
% x0r=Phi'*xn;
% fr=@(xr) Phi'*f(Phi*xr);
% [xDMD1,DMD_time_case1]=implicitEulerSO(Mr,Dr,Br,x0r,fr,del,u,steps);
% yDMD1=Cr*xDMD1./n;
% plot(tSim,yFOM1,'k',tSim,yDMD1,'r--','linewidth',3)

%% Results 
%disp('Case-I results:')
%  %[Error1] = calErrorPS(yFOM1,yROM_pod1,[]);
%  Plot_output_PS(tSim,tSim,yFOM1,yROM_pod1);
%  Plot_error_PS(tSim,tSim,Error1);
%  showTable_PS(Error1,FOM_time_case1,ROM_pod_time1,POD_time1,yFOM1,n,rdefl)