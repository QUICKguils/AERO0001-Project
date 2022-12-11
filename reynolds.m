function re_x = reynolds(U, x, nu)
% REYNOLDS  Calculate the Reynolds number.
%
% Parameters:
%	U: double
%		Free stream velocity [m/s].
%	x: double, optional
%		Characteristic length [m].
%		Default is the wing chord.
%	nu: double, optional
%		Dynamic viscosity of the fluid [Pa*s].
%		Default is the dry air at 15°C.

% Informations about the wind tunnel setup.
% This structure contains the data saved in setup.mat.
lab_set = load('setup.mat');

% Default options.
if nargin < 3
	nu = 15e-6;
end
if nargin < 2
	x  = lab_set.chord;
end

% Compute Reynolds.
re_x = U * x / nu;

end