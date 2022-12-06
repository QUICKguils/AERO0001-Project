function [cl_conv, cd_conv] = hs_convergence(naca_id, cfg, varargin)
% HS_CONVERGENCE  Study the convergence of hess_smith() function.
%
% This function returns the cd and cl calculated by hess_smith() for different
% numbers of used panels. The convergences of cd and cl are then plotted.
%
% Parameters:
%	naca_id: 1x4 char
%		NACA number of the airfoil.
%	cfg: double
%		Type of configuration desired, corresponding to the index of the desired
%		wind tunnel test. The serie of the performed tests is stored in
%		group_5.mat.
%	varargin: 1x1 cell.
%		The cell consist of a char array that holds optional flags 'p' and 'w'.
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file (see write_results()).
%
% Returns:
%	cl_conv, cd_conv: double(1, np_size)
%		Array of lift and drag coefficient computed for different numbers of
%		panels.
%
% This matlab function was written for the project carried out as part of the
% Aerodynamics course (AERO0001-1), academic year 2022-2023.
% author:  Guilain Ernotte <gernotte@student.uliege.be>
% created: 2022-12-03T12:31+02:00

%% Set parameters.

% Unpack the optional flags from varargin.
opts = '';
if ~isempty(varargin)
	opts = varargin{1};
end

% Create sample of panels, cls and cds.
np_sample = 20:20:350;
np_size = numel(np_sample);
cl_conv = zeros(1, np_size);
cd_conv = zeros(1, np_size);

%% Compute lift and drag coefficients.

% Compute cl and cd for each sample item.
for np = 1:np_size
	[~, cd_conv(np), cl_conv(np), ~] = hess_smith(naca_id, np_sample(np), cfg);
end

%% Plot.

if contains(opts, 'p')
	% Plot the convergence of cd and cl.
	figure('WindowStyle', 'docked');
	plot(np_sample, [cl_conv; cd_conv]);
	
	% Dress the figure.
	title("Convergence of the H&S code");
	xlabel("Number of panels");
	ylabel("cl, cd");
	legend("cl", "cd");
	grid;
end

%% Register in external file.

if contains(opts, 'w')
	% Informations about the performed tests in the wind tunnel.
	% This structure contains the data saved in group_5.mat.
	lab_res = load('group_5.mat');

	% Specify the record file name.
	dir = 'Results/';
	filename = strcat( ...
		dir, ...
		'hs-cl_cd_convergence',...
		'-a', num2str(floor(lab_res.AoA(cfg))), ...
		'-v', num2str(floor(lab_res.Uinf(cfg))), ...
		'.csv');
	
	% Gather the data to store.
	ext_conv = [np_sample; cl_conv; cd_conv]';
	
	% Write in external file.
	writematrix(ext_conv, filename);
end
end