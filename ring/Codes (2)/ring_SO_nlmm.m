
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
Jf=@(x) -cos(x) - 100*(cos(x-circshift(x,1)) + cos(x-circshift(x,-1)));
%% Case1: All nodes starting from under-perturbation
x0=0.8*ones(n,1); % initial conditions
steps=k;
disp('Simulating Case-I: All nodes starting from under perturbation..')

[xFOM1,FOM_time_case1]=implicitEulerSO(M,D,B,x0,f,del,u,steps);

yFOM1=C*xFOM1./n;  
%% NLMM SO
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
%%
Phi=Vnlmm1;
Mr=Phi'*M*Phi;
Dr=Phi'*D*Phi;
Br=Phi'*B;
Cr=C*Phi;
x0=Vnlmm1'*x0;
fr=@(xr) Phi'*f(Phi*xr);
[xnlmm1,nlmm_time_case1]=implicitEulerSO(Mr,Dr,Br,x0,fr,del,u,steps);
ynlmm1=Cr*xnlmm1./n;

%% POD
tic
x0=0.8*ones(n,1);
[Vpod,SS,VV]=svd(xFOM1,'econ');
Vpod1=Vpod(:,1:rdefl);
POD_time1=toc;
[xROM_pod,ROM_pod_time1]=implicitEuler_SO_ROM(M,D,B,x0,f,del,u,steps,Vpod1);
yROM_pod1=C*xROM_pod./n;
%% Results 
disp('Case-I results:')
 [Error2] = calErrorPS(yFOM1,yROM_pod1,ynlmm1,[]);
 Plot_output_PS(tSim,tSim,tSim,yFOM1,yROM_pod1,ynlmm1);
 Plot_error_PS(tSim,tSim,Error2);
 