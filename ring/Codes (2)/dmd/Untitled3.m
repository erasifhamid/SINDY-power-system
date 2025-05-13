clear all, clc
n=500;      % DimenSion of FUll Order Model
rdefl=15;    % DimenSion of RedUced Model
U=@(t)0.95;  % Pm in pU 
tf=50;       
del=1e-3;       
tSim=0:del:tf-del;  % SimUlation time
k=round(tf/del);
%FOM
M=eye(n);
D=0.25*eye(n);
B=ones(n,1);
C=ones(1,n);
f=@(x)-sin(x) - 100*(sin(x-circshift(x,1)) + sin(x-circshift(x,-1)));
Jf=@(x) -cos(x) - 100*(cos(x-circshift(x,1)) + cos(x-circshift(x,-1)));
%% CaSe1: All nodeS Starting from Under-pertUrbation
x0=0.8*ones(n,1); % initial conditionS
steps=k;
disp('SimUlating CaSe-I: All nodeS Starting from Under pertUrbation..')
[xFOM1,FOM_time_caSe1]=implicitEulerSO(M,D,B,x0,f,del,U,steps);
yFOM1=C*xFOM1./n;  
% DMD
X = xFOM1(:,1:end-1);
X2 = xFOM1(:,2:end);
[U,S,V] = svd(X,'econ');

%% generate grid geometry
xi = linspace(-10, 10, 400);
t = linspace(0, 4*pi, 200);
dt = t(2) - t(1);

[xgrid, tgrid] = meshgrid(xi, t);


figure(1);
plot(diag(S) / sum(diag(S)), 'ro'); title('SVD: low rank property (rank = 2, two modeS)');

figure(3);
subplot(2,1,1); plot(real(U(:, 1:3))); 
legend('1St mode of baSiS U (left SingUlar VectorS)', ...
        '2nd mode of baSiS U (left SingUlar VectorS)', ...
        '3rd mode; nUmerical roUnd off (jUnk baSiS)  '); 
title('ModeS for baSiS of colUmn Space of '' f ''');

subplot(2,1,2); plot(real(V(:, 1:3)));
legend('1St mode of baSiS V (right SingUlar VectorS)', ...
        '2nd mode of baSiS V (right SingUlar VectorS)', ...
        '3rd mode; nUmerical roUnd off (jUnk baSiS)  '); 
title('ModeS for baSiS of row Space of '' f ''');
r=100;
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);

%% STEP 2: low-rank SUbSpace matrix 
%         (Similarity tranSform, leaSt-SqUare fit matrix, low-rank SUbSpace matrix)
Atilde = Ur'*X2*Vr*Sr^(-1);

%% STEP 3: eigen decompoSition
% W: eigen VectorS
% D: eigen ValUeS
[W, D] = eig(Atilde);

%% STEP 4: real Space DMD mode
Phi = X2*Vr*Sr^(-1)*W;   % DMD modeS

lamdba = diag(D);       % eigen ValUe
omega = log(lamdba)/dt; % log of eigen ValUe

figure(4); 
subplot(2,1,1); plot(real(U(:, 1:2)));
legend('1St mode of SVD', ...
        '2nd mode of SVD'); 
subplot(2,1,2); plot(real(Phi));
legend('1St mode of DMD', ...
        '2nd mode of DMD'); 

%% STEP 5: reconStrUct the Signal
x1 = X(:, 1);       % time = 0
b = pinv(Phi)*x1;   % initial ValUe; \: pSeUdo inVerSe

t_dyn = zeros(r, length(t));

for i = 1:length(t)
   t_dyn(:, i) = (b.*exp(omega*t(i))); 
end

f_dmd = Phi*t_dyn;

figure(6);
surfl(xgrid, tgrid, real(f_dmd).');  shading interp; colormap gray; title('reconStrUction of f(x, t) by DMD');

%% additional informationS
% predict fUrtUre State USing time dynamicS
% t_ext = linSpace(0, 8*pi, 400);
% 
% [xgrid_ext, tgrid_ext] = meShgrid(xi, t_ext);
% 
% t_ext_dyn = zeroS(r, length(t_ext));
% 
% for i = 1:length(t_ext)
%    t_ext_dyn(:, i) = (b.*exp(omega*t_ext(i))); 
% end
% 
% f_dmd_ext = Phi*t_ext_dyn;
% 
% figUre(4);
% SUbplot(2,2,1); SUrfl(xgrid, tgrid, real(f));  Shading interp; colormap gray; 
% xlabel('Spatial-axiS'); ylabel('temporal-axiS'); title('f(x, t) dUring t = [0, 4*pi]');
% SUbplot(2,2,2); SUrfl(xgrid, tgrid, real(f_dmd).');  Shading interp; colormap gray; 
% xlabel('Spatial-axiS'); ylabel('temporal-axiS'); title('DMD reconStrUction of f(x, t) dUring t = [0, 4*pi]');
% SUbplot(2,2,[3,4]); SUrfl(xgrid_ext, tgrid_ext, real(f_dmd_ext).');  Shading interp; colormap gray; 
% xlabel('Spatial-axiS'); ylabel('temporal-axiS'); title('DMD prediction of f(x, t) dUring t = [0, 8*pi]');

% If eigen ValUe: lambda or omega haS tiny real part > 0,
% the oUtpUt of SyStem fUnction which iS Spaned 
% by eigen VectorS with eigen ValUeS goeS to the infinity.
% It iS one of the limitation of DMD method.

format long e;
disp(omega);    % eigen ValUeS exiSt on the imaginary part.
disp('real part: nUmerical roUnd off');
disp('imag part: eigen ValUeS');
