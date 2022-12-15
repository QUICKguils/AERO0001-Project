function out = reynolds(in, x, nu, inverse)
% REYNOLDS  Calculate the Reynolds number or the free stream velocity.
%
% Parameters:
%	in: double
%		Reynolds number of free stream velocity [m/s].
%	x: double, optional
%		Characteristic length [m].
%		Default is the wing chord.
%	nu: double, optional
%		Dynamic viscosity of the fluid [Pa*s].
%		Default is the dry air at 15°C.
%	inverse: boolean, optional
%		Flag to signal that we provide the Reynolds number as input.

% Informations about the wind tunnel setup.
% This structure contains the data saved in setup.mat.
lab_set = load('setup.mat');

% Default options.
if nargin < 2
	x  = lab_set.chord;
end
if nargin < 3
	nu = 15e-6;
end
if nargin < 4
	inverse = false;
end

% Return output.
if ~inverse
	% Compute Reynolds.
	out = in * x / nu;
else
	% Compute free stream velocity.
	out = in * nu / x;
end
end