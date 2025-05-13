% Power-grid synchronization models
% version 1 22-Jan-2015
%
%   Computing parameters for network dynamics models:
%
%   EN_model - Effective network model
%   SP_model - Structure-preserving model
%   SM_model - Synchronous motor model
%
%   Default methods for estimating dynamic parameters:
%
%   default_x_d - Transient reactances x'_{d,i}
%   default_H   - Inertia constants H_i
%   default_D   - Damping coefficients D_i
%
%   Example systems:
%
%   test_system_3gen  - 3-generator test system
%   test_system_10gen - 10-generator (New England) test system
%   test_system_50gen - 50-generator test system
%   poland_case2383wp - Poland power grid
