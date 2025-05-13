% This file is a part of manuscript "Structure Preserving Nonlinear Reduced Order Modeling Technique for Power Systems"
% Authors: Danish Rafiq and Mohammad Abid Bazaz
% Contact: danish_pha2007@nitsri.ac.in

% Solves the Ring Grid Power System Model for three scenarios (as discussed in the paper) in second-order structure
% with n similar generators via SO-NLMM-DEIM & POD using implicit Euler Method
% M*delta_ddot + D* delta_dot = pm-(b*sin(delta_i)+bint[sin(delta_i - delta_i+1)+sin(delta_i - delta_i-1)])

% Selected References:
% [1] Reduced order modelling for transient simulation of power systems
%     using TPWL (Haris Malik, Borzacchiello,Francisco and Diez)
% [2] Nonlinear Moment Matching for the Simulation-Free Reduction of
%     Structural Systems (Maria et al.)

%%
clear all; close all; clc
n=1000;      % Dimension of Full Order Model
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
%% Case1: All nodes starting from under-perturbation
x0=0.8*ones(n,1); % initial conditions
steps=k;
disp('Simulating Case-I: All nodes starting from under perturbation..')
[xFOM1,FOM_time_case1]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM1=C*xFOM1./n;  
% SO_NLMM
K = n; 
thelp = linspace(0,50,K)'; 
SG=1;
v0 = zeros(K,SG);
dqr = zeros(K,SG);
qr = zeros(K,SG);
r = zeros(K,SG);
ddqr = zeros(K,SG);
amp=10; 
omega=pi; 
ph=50;
% Generate qr, dqr, ddqr
for i=1:length(thelp)
    qr(i,1) = amp*sin(omega*thelp(i)+ph);
    dqr(i,1) = amp*cos(omega*thelp(i)+ph);
    ddqr(i,1) = -amp*sin(omega*thelp(i)+ph);
end
r = qr;
optionsMM.rdefl = rdefl;  %deflation of V_NLMM
optionsMM.RelTol = 1e-3; % relative tolerance for residual error: norm(f(xcurr))
optionsMM.AbsTol = 1e-6;
optionsMM.Display = 'none'; % 'none', 'iter'
optionsMM.real = 1;
optionsMM.solver='fsolve';  %NewtonRaphson, fsolve
nlse=15;
tic
[Vnlmm1, nNLSE,~] = SO_nlmm(M,D,B,f,Jf,qr,dqr,ddqr,r,v0,nlse,optionsMM);
NLMMtime1=toc;
%DEIM
Opts.NoSV=10; %DEIM indices
[U_trunc,P,~,~] = funcSnapBasisDEIM(f,xFOM1,Opts);
pre1 = Vnlmm1'*U_trunc/(P'*U_trunc); %r*m
int1=P'*Vnlmm1; 
DEIM_pre_time=toc;
z0_nlmm=Vnlmm1'*x0;
[xROM_nlmm,ROM_nlmm_time1]=ROM_DEIM(M,D,B,z0_nlmm,f,del,u,steps,Vnlmm1,int1,pre1,P);
yROM_nlmm1=C*(Vnlmm1*xROM_nlmm)./n;
% POD
tic
[Vpod,SS,VV]=svd(xFOM1,'econ');
Vpod1=Vpod(:,1:rdefl);
POD_time1=toc;
[xROM_pod,ROM_pod_time1]=implicitEuler_SO_ROM(M,D,B,x0,f,del,u,steps,Vpod1);
yROM_pod1=C*xROM_pod./n;
% Results 
disp('Case-I results:')
 [Error1] = calErrorPS(yFOM1,yROM_pod1,yROM_nlmm1,[]);
 Plot_output_PS(tSim,tSim,tSim,yFOM1,yROM_pod1,yROM_nlmm1);
 Plot_error_PS(tSim,tSim,Error1);
 showTable_PS(Error1,FOM_time_case1,ROM_pod_time1,ROM_nlmm_time1,POD_time1,NLMMtime1,yFOM1,n,rdefl,Opts.NoSV)
%% Case-II: All nodes starting from over-perturbation
x0=1.15*ones(n,1); %(initial conditions) 0, 1, 1.15
steps=k;
disp('Simulating Case-II: All nodes starting from over-perturbation..')
[xFOM2,FOM_time_case2]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM2=C*xFOM2./n;  
% SO_NLMM
K = n; 
thelp = linspace(0,50,K)'; %(0,50,K)
SG=1;
v0 = zeros(K,SG);
dqr = zeros(K,SG);
qr = zeros(K,SG);
r = zeros(K,SG);
ddqr = zeros(K,SG);
amp=10; %10
omega=pi; %pi
ph=50; %50
% Generate qr, dqr, ddqr
for i=1:length(thelp)
    qr(i,1) = amp*sin(omega*thelp(i)+ph);
    dqr(i,1) = amp*cos(omega*thelp(i)+ph);
    ddqr(i,1) = -amp*sin(omega*thelp(i)+ph);
end
r = qr;
optionsMM.rdefl = rdefl;  %deflation of V_NLMM
optionsMM.RelTol = 1e-3; % relative tolerance for residual error: norm(f(xcurr))
optionsMM.AbsTol = 1e-6;
optionsMM.Display = 'none'; % 'none', 'iter'
optionsMM.real = 1;
optionsMM.solver='fsolve';  %NewtonRaphson, fsolve
nlse=15;
tic
[Vnlmm2, nNLSE,~] = SO_nlmm(M,D,B,f,Jf,qr,dqr,ddqr,r,v0,nlse,optionsMM);
NLMMtime2=toc;
%DEIM
Opts.NoSV=10; %10
[U_trunc,P,~,~] = funcSnapBasisDEIM(f,xFOM2,Opts);
pre2 = Vnlmm2'*U_trunc/(P'*U_trunc); %r*m
int2=P'*Vnlmm2; % m*r
DEIM_pre_time=toc;
F_deim=pre2*P';
z0_nlmm=Vnlmm2'*x0;
[xROM_nlmm2,ROM_nlmm_time2]=ROM_DEIM(M,D,B,z0_nlmm,f,del,u,steps,Vnlmm2,int2,pre2,P);
yROM_nlmm2=C*(Vnlmm2*xROM_nlmm2)./n;
% POD
tic
[Vpod,SS,VV]=svd(xFOM2,'econ');
Vpod2=Vpod(:,1:rdefl);
POD_time2=toc;
[xROM_pod2,ROM_pod_time2]=implicitEuler_SO_ROM(M,D,B,x0,f,del,u,steps,Vpod2);
yROM_pod2=C*xROM_pod2./n;
% Results 
disp('Case-II results:')
 [Error2] = calErrorPS(yFOM2,yROM_pod2,yROM_nlmm2,[]);
 Plot_output_PS(tSim,tSim,tSim,yFOM2,yROM_pod2,yROM_nlmm2);
 Plot_error_PS(tSim,tSim,Error2);
 showTable_PS(Error2,FOM_time_case2,ROM_pod_time2,ROM_nlmm_time2,POD_time2,NLMMtime2,yFOM2,n,rdefl,Opts.NoSV)
%% Case-II: All nodes starting from sychronously equillibrumcondition
x0=ones(n,1); %(initial conditions) 0, 1, 1.15
steps=k;
disp('Simulating Case-III: All nodes starting from synchronously equillibrum condition...')
[xFOM3,FOM_time_case3]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM3=C*xFOM3./n;  
% SO_NLMM
K = n; 
thelp = linspace(0,50,K)'; %(0,50,K)
SG=1;
v0 = zeros(K,SG);
dqr = zeros(K,SG);
qr = zeros(K,SG);
r = zeros(K,SG);
ddqr = zeros(K,SG);
amp=10; %10
omega=pi; %pi
ph=50; %50
% Generate qr, dqr, ddqr
for i=1:length(thelp)
    qr(i,1) = amp*sin(omega*thelp(i)+ph);
    dqr(i,1) = amp*cos(omega*thelp(i)+ph);
    ddqr(i,1) = -amp*sin(omega*thelp(i)+ph);
end
r = qr;
optionsMM.rdefl = rdefl;  %deflation of V_NLMM
optionsMM.RelTol = 1e-3; % relative tolerance for residual error: norm(f(xcurr))
optionsMM.AbsTol = 1e-6;
optionsMM.Display = 'none'; % 'none', 'iter'
optionsMM.real = 1;
optionsMM.solver='fsolve';  %NewtonRaphson, fsolve
nlse=15;
tic
[Vnlmm3, nNLSE,~] = SO_nlmm(M,D,B,f,Jf,qr,dqr,ddqr,r,v0,nlse,optionsMM);
NLMMtime3=toc;
%DEIM
Opts.NoSV=10; %10
[U_trunc,P,~,~] = funcSnapBasisDEIM(f,xFOM3,Opts);
pre3 = Vnlmm3'*U_trunc/(P'*U_trunc); %r*m
int3=P'*Vnlmm3; % m*r
DEIM_pre_time=toc;
F_deim=pre3*P';
z0_nlmm=Vnlmm3'*x0;
[xROM_nlmm3,ROM_nlmm_time3]=ROM_DEIM(M,D,B,z0_nlmm,f,del,u,steps,Vnlmm3,int3,pre3,P);
yROM_nlmm3=C*(Vnlmm3*xROM_nlmm3)./n;
% POD
tic
[Vpod,SS,VV]=svd(xFOM3,'econ');
Vpod3=Vpod(:,1:rdefl);
POD_time3=toc;
[xROM_pod3,ROM_pod_time3]=implicitEuler_SO_ROM(M,D,B,x0,f,del,u,steps,Vpod3);
yROM_pod3=C*xROM_pod3./n;
% Results 
disp('Case-III results:')
 [Error3] = calErrorPS(yFOM3,yROM_pod3,yROM_nlmm3,[]);
 Plot_output_PS(tSim,tSim,tSim,yFOM3,yROM_pod3,yROM_nlmm3);
 Plot_error_PS(tSim,tSim,Error3);
 showTable_PS(Error3,FOM_time_case3,ROM_pod_time3,ROM_nlmm_time3,POD_time3,NLMMtime3,yFOM3,n,rdefl,Opts.NoSV)
disp('simulation finished')

