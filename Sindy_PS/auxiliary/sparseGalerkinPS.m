function dy = sparseGalerkinPS(y,ahat)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

yPool = poolDataPS(y',length(y));
dy = (yPool*ahat)';