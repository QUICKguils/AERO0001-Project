function [cd_conv] = hs_convergence(naca_id, opts)
% HS_CONVERGENCE  Study the convergence of hess_smith() function.
%
% This function returns the cd calculated by hess_smith() for different
% numbers of used panels. The convergences of cd are then plotted.
%
% Parameters:
%	naca_id: 1x4 char
%		NACA number of the airfoil.
%	opts: char {'p'|'w'}, optional
%		Optional flags:
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.
%
% Returns:
%	cd_conv: double(3, np_size)
%		Array of lift and drag coefficient computed for different numbers of
%		panels.
%
% This matlab function was written for the project carried out as part of the
% Aerodynamics course (AERO0001-1), academic year 2022-2023.
% author:  Guilain Ernotte <gernotte@student.uliege.be>
% created: 2022-12-03T12:31+02:00

%% Set parameters.

% Set default opts to 'p' (plot, but do not write).
if nargin < 2
	opts = 'p';
end

% Index of the selected w.t. experiment configurations.
cfgs = 1:3;
ncfg = numel(cfgs);

% Create sample of panels, cls and cds.
np_sample = 10:4:600;
np_size = numel(np_sample);
cd_conv = zeros(ncfg, np_size);

%% Compute lift and drag coefficients.

% Compute cd from hess_smith, for 
for cfg = 1:ncfg
	for np = 1:np_size
		[~, cd_conv(cfg, np), ~] = hess_smith(naca_id, np_sample(np), cfg);
	end
end

%% Plot.

if contains(opts, 'p')
	% Plot the convergence of cd and cl.
	figure('WindowStyle', 'docked');
	semilogy(np_sample, cd_conv);
	
	% Dress the figure.
	title("Convergence of the H&S code");
	xlabel("Number of panels");
	ylabel("cd");
	legend("aoa = 5°", "aoa = 10°", "aoa = 15°");
	grid;
end

%% Register in external file.

if contains(opts, 'w')
	% Specify the record file name.
	filename = 'Results/hs-cd_convergence-v16.csv';
	
	% Gather the data to store.
	ext_conv = [np_sample; cd_conv]';
	
	% Write in external file.
	writematrix(ext_conv, filename);
end
end