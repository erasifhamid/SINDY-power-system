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
Jf=@(x) -cos(x) - 100*(cos(x-circshift(x,1)) + cos(x-circshift(x,-1)));
%% Case1: All nodes starting from under-perturbation
x0=0.8*ones(n,1); % initial conditions
steps=k;
disp('Simulating Case-I: All nodes starting from under perturbation..')
[xFOM1,FOM_time_case1]=implicitEulerSO(M,D,B,x0,f,del,u,steps);
yFOM1=C*xFOM1./n;  
% DMD
X = xFOM1(:,1:end-1);
X2 = xFOM1(:,2:end);
[U,S,V] = svd(X,'econ');
%%  Compute DMD (Phi are eigenvectors)
r = 100;  % truncate at 15 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W; 
figure(1)
plot(real(Phi));
legend('1st mode of DMD', ...
        '2nd mode of DMD'); 
[Phi,~] = GramSchmidt(Phi);
%% ROM using DMD
Mr=Phi'*M*Phi;
Dr=Phi'*D*Phi;
Br=Phi'*B;
Cr=C*Phi;
x0r=pinv(Phi)*x0;
fr=@(xr) Phi'*f(Phi*xr);
[xDMD1,DMD_time_case1]=implicitEulerSO(Mr,Dr,Br,x0r,fr,del,u,steps);
yDMD1=Cr*xDMD1./n;
plot(tSim,yFOM1,'k',tSim,yDMD1,'r--','linewidth',3)