function dy = SparseGalerkinie(t,y,tspan,ahat,polyorder,delta,gamma,K)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

f=interp1(tspan,delta,t);

yPool = Pooldataout(y',polyorder,f,gamma,K);
dy = (yPool*ahat)';