function [cl, cd] = wind_tunnel(cfg, opts)
% WIND_TUNNEL  Gather and plot lab results.
%
% The code below was written for the project carried out as part of the
% aerodynamics course (AERO0001-1).
%
% Parameters:
%	cfg: double
%		Type of configuration desired. Corresponds to the index of the
%		desired wind tunnel test. The serie of the performed tests is
%		stored in group_5.mat.
%	opts: char {'p'|'w'}, optional
%		Optional flags:
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.

% Set default opts to an empty char (no plot, no write).
if nargin < 2
	opts = '';
end

% Import the wind tunnel experiment setup.
lab_set = load("setup.mat");
lab_res = load("group_5.mat");
% Unpack.
aoa   = lab_res.AoA(cfg);
v_inf = lab_res.Uinf(cfg);
rho   = lab_set.rho;
c     = lab_set.chord;

% % Wind tunnel corrections.
% K1 = 0.52; % TODO: make sure of this coef.
% wt_width = 2.5;
% wt_height = 1.8;
% wt_S = wt_width * wt_height;
% wg_L = 1; % Approx of the wing length [m].
% wg_t = lab_set.chord * 0.18; % Max tickness of the airfoil [m].
% wg_V = lab_set.chord * wg_t * wg_L;
% eps = K1 * wg_V / wt_S^(3/2);
% % -> 2e-3 -> correction factor is negligible.

%% Compute, plot and write in external file the cp vs chord.

xc = lab_set.coord_taps(1, :)/lab_set.chord;
cp = lab_res.p / (0.5 * lab_set.rho * lab_res.Uinf(cfg)^2);

if contains(opts, 'p')
	plot_cp(xc, cp);
end
if contains(opts, 'w')
	write_results(xc, cp);
end

%% Computation of lift and drag coefficients.

% Indexes of the pressure taps that faces each other on the upper
% and lower side of the airfoil.
taps_up  = flip([1:12, 14:19]);
taps_low = [21:33, 35:39];

% Pressures on the upper and lower sides.
p_up = lab_res.p(cfg, taps_up);
p_low = lab_res.p(cfg, taps_low);

% x-coordinate of the taps. They are the same on upper and lower side.
x_taps = lab_set.coord_taps(1, taps_up);  % or taps_low. Equivalent.

% Differences between adjacent x-coordinates.
dx = diff([x_taps,  c]);
% Differences between adjacent y-coordinates.
dy_up  = diff([lab_set.coord_taps(2, taps_up),  0]);
dy_low = diff([lab_set.coord_taps(2, taps_low), 0]);  % It is just the opposite.

% Axial and normal forces.
Nprime = trapz(x_taps, p_low-p_up);
Aprime = trapz(x_taps, p_up.*(dy_up./dx) - p_low.*(dy_low./dx));

% Lift and drag forces.
Lprime = Nprime*cosd(aoa) - Aprime*sind(aoa);
Dprime = Nprime*sind(aoa) + Aprime*cosd(aoa);

% Lift and drag coefficients.
cl = Lprime / (0.5*rho*v_inf^2*c);
cd = Dprime / (0.5*rho*v_inf^2*c);

%% 5. Local function definitions.

	function plot_cp(xc, cp)
		% PLOT_CP  Plot the chordwise distribution of the cp.

		figure('WindowStyle', 'docked');
		hold on;

		% Upper side.
		plot( ...
			xc(1:floor(end/2)),      ...
			cp(cfg, 1:floor(end/2)), ...
			'color', 'red',          ...
			'Marker', 'x',           ...
			'Linewidth', 1);

		% Lower side.
		plot( ...
			xc(ceil(end/2):end),      ...
			cp(cfg, ceil(end/2):end), ...
			'color', 'blue',          ...
			'Marker', 'x',            ...
			'Linewidth', 1);

		% Dress the plot.
		title('Distribution of the pressure coefficient Cp along the chord');
		xlabel('x/c');
		ylabel('cp');
		grid;
		legend('Upper surface', 'Lower surface');
		set(gca, 'YDir', 'reverse')
	end

	function write_results(xc, cp)
		% WRITE_RESULTS  Write desired results in external files.

		% Specify the record file name.
		filename = strcat( ...
			'Results/', ...
			'lab-cp', ...
			'-a', num2str(floor(aoa)), ...
			'-v', num2str(floor(v_inf)));

		% Gather the data to store.
		plot_up = [ ...
			xc(1:floor(end/2)); ...
			cp(config, 1:floor(end/2))]';
		plot_low = [ ...
			xc(ceil(end/2):end); ...
			cp(config, ceil(end/2):end)]';

		% Write in external file.
		writematrix(plot_up,  strcat(filename, '-up.csv'));
		writematrix(plot_low, strcat(filename, '-low.csv'));
		writematrix(lab_set.coord_taps', 'Results/pressure_taps.csv');
	end
end
