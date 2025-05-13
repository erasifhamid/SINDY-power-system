function mpc = poland_case2383wp
%POLAND_CASE2383WP  Polish system (winter 1999-2000 peak).
%   mpc = POLAND_CASE2383WP just loads the power flow data case2383wp
%   provided in the MATPOWER package and sets the reference frequency to be
%   50 Hz. When feeding mpc to the functions EN_model, SP_model, or
%   SM_model, the dynamic parameters will be estimated according to the
%   method used in Ref. [1,2].
%
%   References:
%
%   [1] T. Nishikawa and A. E. Motter, Comparative analysis of existing
%   models for power-grid synchronization, New J. Phys. 17, 015012 (2015).
%
%   [2] A. E. Motter, S. A. Myers, M. Anghel, and T. Nishikawa, Spontaneous
%   synchrony in power-grid networks, Nat. Phys. 9, 191-197 (2013).
%
%   See also default_D, default_H, default_x_d, EN_model, SP_model,
%   SM_model
%
%   Comments from the original dataset 'case2383wp' in MATPOWER:
%
%   This case represents the Polish 400, 220 and 110 kV networks during
%   winter 1999-2000 peak conditions. It is part of the 7500+ bus Europen
%   UCTE system. To decrease the number of buses, the tie lines to foreign
%   networks were replaced by artificial load or generator buses (180-186).
%   Multiple generators at a bus have been aggregated. Generators that are
%   not centrally dispatchable in the Polish energy market are given a cost
%   of zero.
%  
%   This data was graciously provided by, and is distributed with the
%   permission of, Roman Korab <roman.korab@polsl.pl>.

%
% Copyright (C) 2015  Takashi Nishikawa
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at
% your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
% USA.

%   Last modified by Takashi Nishikawa on 1/27/2015

mpc = loadcase('case2383wp.m');
mpc.ref_freq = 50;
