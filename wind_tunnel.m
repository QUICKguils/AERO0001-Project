function [cl, cd] = wind_tunnel(cfg, opts)
% WIND_TUNNEL  Gather and plot lab results.
%
% The code below was written for the project carried out as part of the
% aerodynamics course (AERO0001-1).

% TODO:
% - wind tunnel corrections. Apparently, the corrections are quite negligible.
% Caculate the eps factor, and say in the report that we would not take into
% account this negligible correction.

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

% Reynolds number
% Rough: v16 -> 4.81e5
%        v25 -> 7.44e5

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

% Pressures on the upper and lower sides.
p_up = flip(lab_res.p(cfg, 1:19));
p_low = lab_res.p(cfg, 21:end);

dx_up = diff(flip([c, lab_set.coord_taps(1, 1:19)]));
dy_up = diff(flip([0, lab_set.coord_taps(2, 1:19)]));
dx_low = diff([lab_set.coord_taps(1, 21:end), c]);
dy_low = diff([lab_set.coord_taps(2, 21:end), 0]);

% Axial and normal forces.
Nprime = trapz(lab_set.coord_taps(1, 21:end), p_low-p_up);
Aprime = trapz(lab_set.coord_taps(1, 21:end), p_up.*(dy_up./dx_up) - p_low.*(dy_low./dx_low));

% Lift and drag forces.
Lprime = Nprime*cosd(aoa) - Aprime*sind(aoa);
Dprime = Nprime*sind(aoa) + Aprime*cosd(aoa);

% Lift and drag coefficients.
cl = Lprime / (0.5*rho*v_inf^2*c);
cd = Dprime / (0.5*rho*v_inf^2*c);

%% 5. Local function definitions.

	function plot_cp(xc, cp)
		% Pressure coefficient graph.
		figure('WindowStyle', 'docked'); hold on; grid;
		plot( ...
			xc(1:floor(end/2)),      ...
			cp(cfg, 1:floor(end/2)), ...
			'color', 'red',          ...
			'Marker', 'x',           ...
			'Linewidth', 1);
		plot( ...
			xc(ceil(end/2):end),      ...
			cp(cfg, ceil(end/2):end), ...
			'color', 'blue',          ...
			'Marker', 'x',            ...
			'Linewidth', 1);
		xlabel('x/c');
		ylabel('Cp');
		legend('Upper surface', 'Lower surface');
		title('Distribution of the pressure coefficient Cp along the chord');
		set(gca, 'YDir', 'reverse')
	end

	function write_results(xc, cp)
		% Uncomment this section to write plotting data into external file.
		plot_up = [ ...
			xc(1:floor(end/2)); ...
			cp(config, 1:floor(end/2))]';
		plot_low = [ ...
			xc(ceil(end/2):end); ...
			cp(config, ceil(end/2):end)]';
		writematrix(plot_up, 'Results/lab-cp-a15v25up.csv');
		writematrix(plot_low, 'Results/lab-cp-a15v25low.csv');
		writematrix(lab_set.coord_taps', 'Results/pressure_taps.csv');
	end
end