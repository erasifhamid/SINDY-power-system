function [V, nNLSE, v0] = SO_nlmm(M,D,B,f,Jf,qr,dqr,ddqr,r,v0,K,varargin)
% Second-Order NLMM (MATLAB version)
%Syntax: 
%       V = SO_NLMM(f,s,r,xr,v0)
%       V = SO_NLMM(f,s,r,xr,v0,Opts)
% 
% Description:
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       -f:     right-hand side of nonlinear system ODE f(x,u)
%       -s:     row vector of time-discrete values of shift function s(xr)
%               for every moment and time point, values in s, r and xr have
%               to correspond
%       -r:     row vector of time-discrete values of input function r(xr)
%               for every moment and time point, values in s, r and xr have
%               to correspond
%       -xr:    row vector of discrete time values of signals for every
%               moment, corresponding with s and r
%       -v0:    first guess or matrix of first guesses for each time step of
%               every moment for Newton-iteration
%       *Optional Input Arguments:*
%       -Opts:         a structure containing following options
%           -.SpecifyObjectiveGradient:
%           -.MaxIter:
%           -.TolX:
%           -.rdefl:    desired deflated order
%
% Output Arguments:
%       -V:     orthonormal nonlinear moment matching projection basis
%
% Examples:
%
% See Also: 
%       SignalGeneratorTemplate
%
% References:
%       * *[1] Grimme (1997)*, Krylov projection methods for model reduction
%       * *[2] Antoulas (2010)*, Interpolatory model reduction of
%              large-scale dynamical Systems
%
%% parse Options
if nargin >= 7
    Opts = varargin{1};
end

% Options for nlmm
Def.solver = 'fsolve'; %solver: 'NewtonRaphson' or 'fsolve' 
Def.rdefl = [];
Def.real = 0;
Def.Orthogonalization = 'none'; %'GramSchmidt','QR' or 'none' (recommended)

% Options for 'fsolve' AND 'NewtonRaphson' 
Def.MaxIter = 400; 
Def.Display = 'none';

% Options for 'fsolve'
Def.TolX = 1e-6;
Def.TolFun = 1e-6;
Def.Algorithm = []; % 'trust-region-dogleg', 'trust-region-reflective', 'levenberg-marquardt'
Def.SpecifyObjectiveGradient = true; % provides Jacobian for fsolve

if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%%
% dimension of system
n = size(v0,1);

% Number of SGs
nSG = size(v0,2);

% number of time snapshots per signal generator
%K= size(v0,1);  % chose the number of NlSEs each of dimension nhere
%K=28;

% Number of total NLSEs
nNLSE = sum(K);

% preallocation
Vraw = zeros(n,nNLSE);

%% calculate basis
%initiate iSG
iSG = 1;

while iSG <= nSG % for iSG = 1 : nSG
    %if find(xr{iSG}==0) % determine if signal xr(t) contains zeros
    %    nonzero = false;
    %else
    %    nonzero = true;
    %end
    
    for kSnap = 1 : K(iSG)     
         qrik = qr(kSnap);
         vik =  v0(kSnap);
        
        %if nonzero % nonzero = true
        %    fun = @(v) f(v,r{iSG}(:,kSnap)) - v*s{iSG}(kSnap)/xr{iSG}(kSnap);
         %   Jac = @(v) Jf(v,r{iSG}(:,kSnap)) - speye(n)*(s{iSG}(kSnap)/xr{iSG}(kSnap));
        %else % nonzero = false
            %fun = @(v) f(v*xr{iSG}(kSnap),r{iSG}(:,kSnap)) - v*s{iSG}(kSnap);
            %Jac = @(v) Jf(v*xr{iSG}(kSnap),r{iSG}(:,kSnap))*xr{iSG}(kSnap) - speye(n)*s{iSG}(kSnap);
            fun= @(v) M*v * ddqr(kSnap) + D*v*dqr(kSnap) + f(v*qrik) - B*r(kSnap);
            Jac= @(v) M*ddqr(kSnap) + D*dqr(kSnap) + Jf(v*qrik)*qrik;
        %end
        
        if (strcmp(Opts.solver,'fsolve'))
            Vraw(:,(iSG-1)*K(iSG)+kSnap) = fsolve(inputfsolve(fun,Jac),v0(kSnap),Opts); % changed from v0(kSnap)
        else 
            Vraw(:,(iSG-1)*K(iSG)+kSnap) = NewtonRaphson(fun,v0(kSnap),Jac,Opts);
        end
        
        %%-- Real subspace: split complex conjugate columns into real and imag
        if ~isreal(v0(iSG))
            Vraw(:,(iSG*K(iSG)+kSnap)) = conj(Vraw(:,((iSG-1)*K(iSG)+kSnap)));
            if Opts.real
                Vraw(:,(iSG*K(iSG)+kSnap)) = imag(Vraw(:,((iSG-1)*K(iSG)+kSnap)));
                Vraw(:,((iSG-1)*K(iSG)+kSnap)) = real(Vraw(:,((iSG-1)*K(iSG)+kSnap)));
            end 
        end
        
        %%-- Orthogonalization
        if strcmp(Opts.Orthogonalization,'GramSchmidt')
            Vraw = gramSchmidt((iSG-1)*K(iSG)+kSnap, Vraw);
        end
        if strcmp(Opts.Orthogonalization,'QR')
            Vraw = qr(Vraw(:,((iSG-1)*K(iSG)+kSnap)));
        end
        
        if kSnap < K(iSG)
            v0 = Vraw(:,(iSG-1)*K(iSG)+kSnap);
            if ~isreal(v0(iSG))
                v0(kSnap+1) = Vraw(:,(iSG*K(iSG)+kSnap));
            end
        end
    end
    
    if ~isreal(v0(iSG))
        iSG = iSG + 2; % due to complex conjugate pair, we can jump one SG
    else
        iSG = iSG + 1;
    end
    
end

disp(['highest useful reduced order is ' num2str(rank(Vraw))])
% if ~isempty(Opts.rdefl)
%     if Opts.rdefl>rank(Vraw)
%        Opts.rdefl=rank(Vraw);
%        warning(['Setting rdefl to the highest useful reduced order ' num2str(rank(Vraw))])
%     end
% end

%%-- Deflation
%set reduced order via singular value decomposition
if ~isempty(Opts.rdefl) 
    [V] = basisRed(Vraw,Opts.rdefl);
else
    [V] = basisRed(Vraw);
end

end

%%-------------------------- AUXILIARY FUNCTIONS --------------------------
function [fun,Jac] = inputfsolve(fun,Jac)
end

function [V, TRv] = gramSchmidt(jCol, V, GramSchmOpts)
% [V, TRv, W, TLw] = gramSchmidt(jCol, V, W, Opts)

%   Gram-Schmidt orthonormalization
%   Input:  jCol:  Column to be treated
%           V, W:  Krylov-Subspaces
%   Output: V, W:  orthonormal basis of Krylov-Subspaces
%           TRv, TLw: Transformation matrices

% Opts struct
Def.orth = '2mgs';
Def.dgksTol = 1e-12;

if ~exist('GramSchmOpts','var') || isempty(GramSchmOpts)
    GramSchmOpts = Def;
else
    GramSchmOpts = parseOpts(GramSchmOpts,Def);
end

IP = @(x,y) (x.'*y);

TRv=eye(size(V,2));
% TLw=eye(size(V,2));
if jCol>1
    switch GramSchmOpts.orth
        case 'dgks'
            % iterates standard gram-schmidt
            orthError=1;
            count=0;
            while(orthError>GramSchmOpts.dgksTol)
                h=IP(V(:,1:jCol-1),V(:,jCol));
                V(:,jCol)=V(:,jCol)-V(:,1:jCol-1)*h;
                TRv(:,jCol)=TRv(:,jCol)-TRv(:,1:jCol-1)*h;
%                 if hermite
%                     h=IP(W(:,1:jCol-1),W(:,jCol));
%                     W(:,jCol)=W(:,jCol)-W(:,1:jCol-1)*h;
%                     TLw(:,jCol)=TLw(:,jCol)-TLw(:,1:jCol-1)*h;
%                 end
                orthError=norm(IP([V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))],...
                    [V(:,1:jCol-1),V(:,jCol)/sqrt(IP(V(:,jCol),V(:,jCol)))])-speye(jCol),'fro');
                if count>50 % if dgksTol is too small, Matlab can get caught in the while-loop
                    error('Orthogonalization of the Krylov basis failed due to the given accuracy.');
                end
                count=count+1;
            end
        case 'mgs'
            for iCol=1:jCol-1
              h=IP(V(:,jCol),V(:,iCol));
              V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
              TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
%               if hermite
%                 h=IP(W(:,jCol),W(:,iCol));
%                 W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
%                 TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
%               end 
            end
        case '2mgs'
            for k=0:1
                for iCol=1:jCol-1
                  h=IP(V(:,jCol),V(:,iCol));
                  V(:,jCol)=V(:,jCol)-V(:,iCol)*h;
                  TRv(:,jCol)=TRv(:,jCol)-h*TRv(:,iCol);
%                   if hermite
%                     h=IP(W(:,jCol),W(:,iCol));
%                     W(:,jCol)=W(:,jCol)-W(:,iCol)*h;
%                     TLw(:,jCol)=TLw(:,jCol)-h*TLw(:,iCol);
%                   end 
                end
            end
        otherwise
            error('Opts.orth is invalid.');
    end  
end

% normalize new basis vectorm
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
h = sqrt(IP(V(:,jCol),V(:,jCol)));
V(:,jCol)=V(:,jCol)/h;
TRv(:,jCol) = TRv(:,jCol)/h;
% if hermite
%     h = sqrt(IP(W(:,jCol),W(:,jCol)));
%     W(:,jCol)=W(:,jCol)/h;
%     TLw(:,jCol) = TLw(:,jCol)/h;
% end
end

