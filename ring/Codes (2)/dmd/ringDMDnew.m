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
n=500;      % Dimension of Full Order Model
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
%Jf=@(x) -cos(x) - 100*(cos(x-circshift(x,1)) + cos(x-circshift(x,-1)));
% Case1: All nodes starting from under-perturbation
x0=1.15*ones(n,1); % initial conditions
steps=k;
disp('Simulating Case-I: All nodes starting from under perturbation..')
[xFOM1,FOM_time_case1]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM1=C*xFOM1./n;
figure
plot(tSim,yFOM1,'linewidth',2)
%% DMD
tic
X1 = xFOM1(:, 1:end-1);
X2 = xFOM1(:, 2:end);
[U, S, V] = svd(X1, 'econ');
%  Compute DMD (Phi are eigenvectors)
r = 16;      % rank-r truncation
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);
Atilde = Ur'*X2*Vr*Sr^(-1);
[W,eigs] = eig(Atilde);
[W,~] = mgson(W);
Phi=Ur*W;
DMD_time1=toc;
%Phi = X2*Vr*Sr^(-1)*W;
figure
plot(real(Phi))
legend('DMD Modes')
title('Plot of real values of DMD Modes')
%% ROM using DMD
Mr=Phi'*M*Phi;
Dr=Phi'*D*Phi;
Br=Phi'*B;
Cr=C*Phi;
xr=Phi'*x0;
% x0r=pinv(Phi)*x0;
fr=@(x0r) Phi'*f(Phi*x0r);
[xDMD1,DMD_time_case1]=implicitEulerSO(Mr,Dr,Br,xr,fr,del,u,steps);
yDMD1=Cr*xDMD1./n;
plot(tSim,yFOM1,'k',tSim,yDMD1,'r--','linewidth',3)
%% POD
tic
[Vpod,SS,VV]=svd(xFOM1,'econ');
Vpod1=Vpod(:,1:rdefl);
POD_time1=toc;
[xROM_pod,ROM_pod_time1]=implicitEuler_SO_ROM(M,D,B,x0,f,del,u,steps,Vpod1);
yROM_pod1=C*xROM_pod./n;
%% ERROR
error_states_DMD=norm(xFOM1-Phi*xDMD1)/norm(xFOM1)
error_states_POD=norm(xFOM1-Vpod1*xROM_pod)/norm(xFOM1)
%% plotting
figure
subplot(2,2,1)
plot(tSim,yFOM1,'b','linewidth',3)
legend('FULL ORDER MODEL')
title('FULL ORDER MODEL')
subplot(2,2,2)
plot(tSim,yDMD1,'r','linewidth',3)
legend('Reduced Order model with DMD')
title('Reduced Order model with DMD')
subplot(2,2,3)
plot(tSim,yROM_pod1,'g','linewidth',3)
legend('Reduced Order model with POD')
title('Reduced Order model with POD')
subplot(2,2,4)
plot(abs(yFOM1-yROM_pod1),'linewidth',3)
hold on
plot(abs(yFOM1-yDMD1),'linewidth',3)
legend('Error with POD','Error with DMD')
Plot_output_PS(tSim,tSim,tSim,yFOM1,yROM_pod1,yDMD1);
title('Comparison between POD and DMD')
%%

%% Results 

%disp('Case-I results:')
%  %[Error1] = calErrorPS(yFOM1,yROM_pod1,[]);
%  Plot_output_PS(tSim,tSim,yFOM1,yROM_pod1);
%  Plot_error_PS(tSim,tSim,Error1);
%  showTable_PS(Error1,FOM_time_case1,ROM_pod_time1,POD_time1,yFOM1,n,rdefl)