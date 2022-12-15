% TODO: reformat the file registration. Badly implemented.

function cl_aoa(v, opts)
% CL_VS_AOA  Comparison of cl(aoa) between different methods.
%
% Parameters:
%	v: double {16|25}, optional
%		Free stream velocity [m/s]. Default is 16 m/s.
%	opts: char {'p'|'w'}, optional
%		Optional flags:
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file.

%% Parameter settings.

% Set defaults and assert inputs.
if nargin == 0
	v = 16;
end
if nargin <= 1
	opts = 'p';
end
if nargin >= 1
	assert( ...
		v == 16 || v == 25, ...
		"Wind tunnel tests have only been performed for 16 m/s and 25 m/s.")
end

% Determine the lab configuration from the desired speed v.
if v == 16
	cfg = 1:3;
elseif v == 25
	cfg = 4:6;
end

% Import the wind tunnel experiment setup.
lab_res = load('group_5.mat');

aoas_lab = lab_res.AoA(cfg);
aoas_num = -5:1:20;

% init cl.
cl_hs    = zeros(size(aoas_num));
cl_lab   = zeros(size(aoas_lab));

%% Lift coefficient computations.

% cl hs (our panel code).
for aoa = 1:length(aoas_num)
	[~, ~, cl_hs(aoa), ~] = hess_smith('0018', 200, [aoas_num(aoa), v]);
end

% cl lab.
for aoa = 1:length(aoas_lab)
	[cl_lab(aoa), ~] = wind_tunnel(cfg(aoa));
end
% [cl_lab, ~] = arrayfun(@wind_tunnel, cfg(aoas_lab));

% cl from conformal mapping.
tc = 0.18;
cl_cm = @(x) 2*pi * (1 + 4*tc/(3*sqrt(3))) * sin(deg2rad(x));

%% Results plotting.

% Plot these cl.
if contains(opts, 'p')
	% Instantiate the figure.
	figure('WindowStyle','docked');
	hold on;

	% Plot all the methods.
	plot(aoas_num, cl_hs, 'Marker','x');
	plot(aoas_lab, cl_lab, 'Marker','^');
	fplot(cl_cm, [-5, 20]);

	% Dress the plot.
	title(['Cl vs aoa for v = ', num2str(v), ' m/s']);
	xlabel("aoa (deg)");
	ylabel("cl");
	grid;
	legend( ...
		"Matlab code","Wind tunnel", "Conformal mapping", ...
		'Location', 'northwest');
end

%% Results registration.

if contains(opts, 'w')
	% Specify the record file name.
	filename = strcat( ...
		'Results/', ...
		'lab_hs_the-cl_aoa', ...
		'-v', num2str(v), ...
		'.csv');

    % Open the file for writing.
    f = fopen(filename, 'w');

    % Write the data to the file.
    fprintf(f, 'aoa, cl_matlab, cl_conformal_mapping, cl_wind_tunnel\n');
    for aoa = 1:length(aoas_num)
        fprintf(f, '%d, %.4f, %.4f, %.4f, %.4f, ', ...
			aoas_num(aoa), cl_hs(aoa), cl_cm(aoas_num(aoa)));
		% FIXME: this is not robust at all.
		if any(aoas_num(aoa) == aoas_lab)
			fprintf(f, '%.4f', cl_lab((aoa-6)/5));
		end
		fprintf(f, '\n');
    end

    % Close the file.
    fclose(f);
end

end