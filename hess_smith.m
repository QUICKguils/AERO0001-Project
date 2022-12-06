function [cp, cd, cl, cl_KJ] = hess_smith(naca_id, npanel, cfg, varargin)
% HESS_SMITH  Implementation of the Hess & Smith straight panel method.
%
% 1. Parameter settings.
% 2. Airfoil and panel geometry.
% 3. System assembly and solving.
% 4. Derivation and saving of relevant results.
% 5. Local function definitions.
%
% Parameters:
%	naca_id: 1x4 char
%		NACA number of the airfoil.
%	npanel: double
%		Number of panels used.
%	cfg: double
%		Type of configuration desired, corresponding to the index of the desired
%		wind tunnel test. The serie of the performed tests is stored in
%		group_5.mat.
%	varargin: 1x1 cell.
%		The cell consist of a char array that holds optional flags 'p' and 'w'.
%		'p' -> Enable plots creation.
%		'w' -> Write data in external file (see write_results()).
%
% This matlab function was written for the project carried out as part of the
% Aerodynamics course (AERO0001-1), academic year 2022-2023.
% author:  Guilain Ernotte <gernotte@student.uliege.be>
% created: 2022-10-15T14:48+02:00

%% 1. Parameter settings.

% Unpack the optional flags from varargin.
opts = '';
if ~isempty(varargin)
	opts = varargin{1};
end

% Informations about the wind tunnel setup.
% This structure contains the data saved in setup.mat.
lab_set = load('setup.mat');
% Informations about the performed tests in the wind tunnel.
% This structure contains the data saved in group_5.mat.
lab_res = load('group_5.mat');

% General parameters.
c       = lab_set.chord;      % Chordwise length of the wing [m].
rho     = lab_set.rho;        % Air density [kg/m³].
v_inf   = lab_res.Uinf(cfg);  % Free stream velocity [m/s].
aoa     = lab_res.AoA(cfg);   % Angle of attack [°].
aoa_rad = deg2rad(aoa);

% Assertions on parameters value.
assert(~mod(npanel, 2), ...
	"Number of panels must be even.");
assert(all([c, npanel, v_inf] >= 0), ...
	"Aberrant input value detected (negative quantity).");
assert(length(naca_id) == 4, ...
	"This script only supports NACA 4-digits airfoils.");

%% 2. Airfoil and panel geometry.

% Aifoil and camber line.
[X, Y, xh, Y_bar, dx, dy, l, theta] = airfoil_geometry();

% Panel control points, located @ mid-point of each panel.
x_bar = movmean(X, [0, 1]);
y_bar = movmean(Y, [0, 1]);
% Remove irrelevant last value.
x_bar = x_bar(1:end-1);
y_bar = y_bar(1:end-1);

% Plot the airfoil, camber line and panels.
if contains(opts, 'p')
	plot_airfoil(X, Y, xh, Y_bar, x_bar, y_bar);
end

%% 3. System assembly and solving.

[M, N, D, F, C, G, M_11] = build_system(X, Y, x_bar, y_bar, l, theta);

res = M\N;

% Fictive sources associated to each panel.
q = res(1:end-1);
% Vorticity around the airfoil.
gamma = res(end);

%% 4. Derivation and saving of relevant results.

[cp, cd, cl, cl_KJ] = derive_results( ...
	q, gamma, x_bar, dx, dy, l, theta, D, F, C, G, M_11);

% Plot the cp distribution along the chord.
if contains(opts, 'p')
	plot_cp(x_bar, cp);
end

% Write desired results in external files.
if contains(opts, 'w')
	write_results(X, Y, x_bar, cp, c);
end

%% 5. Local function definitions.

	function [X, Y, xh, Y_bar, dx, dy, l, theta] = airfoil_geometry()
		% AIRFOIL_GEOMETRY  Compute the relevant airfoil geometrical parameter.
		%
		% Returns:
		%	X, Y: double(1, npanel+1)
		%		x and y coordinates of the airfoil surface points.
		%	xh: double(1, npanel/2 + 1)
		%		x coordinate along the x axis.
		%	Y_bar: double(1, npanel/2 + 1)
		%		y coordinate of the camber line.
		%	dx, dy: double(1, npanel)
		%		Differences bewteen, respectively, x and y coordinates of
		%		succesive airfoil surface points.
		%	l: double(1, npanel)
		%		Panels length.
		%	theta: double(1, npanel)
		%		Panels inclination w.r.t. the x-axis.

		% Airfoil parameters.
		eps = str2double(naca_id(1))  /100;  % Maximal camber ratio.
		p   = str2double(naca_id(2))  /10;   % Location of maximal camber from LE.
		tau = str2double(naca_id(3:4))/100;  % Thickness ratio.

		% panel definition: x-component boundary points.
		xi = linspace(0, 2*pi, npanel+1);
		x  = c/2 * (cos(xi) + 1);

		% x-half : x coordinate along the x axis.
		xh = x(npanel/2 +1 : end);

		% NACA definition of thickness.
		T = 10 * tau * c * ( ...
			  0.2969 * sqrt(xh/c)    ...
			- 0.1260 *     (xh/c)    ...
			- 0.3537 *     (xh/c).^2 ...
			+ 0.2843 *     (xh/c).^3 ...
			- 0.1015 *     (xh/c).^4);

		% NACA definition of camber line.
		inf_pc = xh < p*c;
		Y_bar = [ ...
			eps*   xh( inf_pc)  /    p ^2 .* (  - xh( inf_pc)/c + 2*p), ...
			eps*(c-xh(~inf_pc)) / (1-p)^2 .* (1 + xh(~inf_pc)/c - 2*p)];

		% First-order approximation of the camber line derivative.
		dY_bar = diff(Y_bar)./diff(xh);
		% Add the backward difference for the last point (@ TE).
		dY_bar = [dY_bar, dY_bar(end)];
		theta_camber = atan(dY_bar);

		% Coordinates of the lower and upper surface.
		X_u = xh    - T/2 .* sin(theta_camber);
		X_l = xh    + T/2 .* sin(theta_camber);
		Y_u = Y_bar + T/2 .* cos(theta_camber);
		Y_l = Y_bar - T/2 .* cos(theta_camber);

		% Gather upper and lower surface coordinate points. Unlike the
		% x-component boundary points, the (X, Y) coordinate points run
		% clockwise around the airfoil, starting and finishing at the TE. We
		% trim the first value of X_u (i.e. the coordinate point of the LE) to
		% not repeat it, as this point is already taken in account by the last
		% element of X_l.
		X = [flip(X_l), X_u(2:end)];
		Y = [flip(Y_l), Y_u(2:end)];

		% Diffs between successive boundary points.
		dx = diff(X);
		dy = diff(Y);
		% Panels length.
		l = sqrt(dx.^2 + dy.^2);
		% Panels inclination w.r.t. the x-axis.
		theta = atan2(dy, dx);
	end

	function plot_airfoil(X, Y, xh, Y_bar, x_bar, y_bar)
		% PLOT_AIRFOIL  Plot the airfoil, camber line and panels.

		figure('WindowStyle', 'docked');
		hold on;

		% Airfoil.
		plot(X, Y, 'marker', '*');
		% Camber line.
		plot(xh, Y_bar);
		% Panel locations.
		plot(x_bar, y_bar, 'linestyle', 'none', 'marker', 'o');

		% Dress the plot.
		title(['NACA-', naca_id, ' airfoil preview']);
		xlabel('X/m'); ylabel('Y/m');
		axis equal; grid;
		legend( ...
			'Panel boundary points.', ...
			'Camber line.', ...
			'Panel control points.');
	end

	function [M, N, D, F, C, G, M_11] = build_system(X, Y, x_bar, y_bar, l, theta)
		% BUILD_SYSTEM  Build the Hess & Smith linear system.
		%
		% Parameters:
		%
		% Returns:
		%	D, F, C, G, M_11: double(npanel, npanel)
		%		See theoretical course 5 for implementation details.
		%	M, N: double(npanel+1, npanel+1), double(npanel+1, 1)
		%		M and N matrices of the linear system to solve.

		% Matrices preallocation.
		A = zeros(npanel, npanel);
		B = zeros(npanel, npanel);
		C = zeros(npanel, npanel);
		D = zeros(npanel, npanel);
		E = zeros(npanel, npanel);
		F = zeros(npanel, npanel);
		G = zeros(npanel, npanel);

		% Matrices construction.
		for i = 1:npanel
			for j = 1:npanel
				A(i,j) = ...
					-(x_bar(i)-X(j)) * cos(theta(j)) ...
					-(y_bar(i)-Y(j)) * sin(theta(j));
				E(i,j) = ...
					 (x_bar(i)-X(j)) * sin(theta(j)) ...
					-(y_bar(i)-Y(j)) * cos(theta(j));
				C(i,j) = sin(theta(i)-theta(j));
				D(i,j) = cos(theta(i)-theta(j));
				B(i,j) = (x_bar(i)-X(j))^2 + (y_bar(i)-Y(j))^2;
				F(i,j) = log(1 + (l(j)^2 + 2*A(i,j)*l(j)) / B(i,j));
				G(i,j) = atan2((E(i,j)*l(j)), A(i,j)*l(j)+B(i,j));
			end
		end

		% Contruction of matrix M.

		M_11 = C.*F/2 - D.*G;
		M_12 = D.*F/2 + C.*G;

		% Force the values on the M_11 and M_12 diagonals.
		M_11 = M_11 - diag(diag(M_11)) + pi*eye(size(M_11));  % pi on M_11 diagonal
		M_12 = M_12 - diag(diag(M_12));                       % 0  on M_12 diagonal

		M_12 = sum(M_12,2);

		M_21 = ...
			  D(end,:).*F(end,:)/2 + C(end,:).*G(end,:) ...
			+ D(1,  :).*F(1,  :)/2 + C(1,  :).*G(1,  :);
		M_22 = - sum( ...
			  C(end,:).*F(end,:)/2 - D(end,:).*G(end,:) ...
			+ C(1  ,:).*F(1,  :)/2 - D(1,  :).*G(1,  :));

		M = 1/(2*pi) * [M_11, M_12; M_21, M_22];

		% Construction of column vector N.

		N = v_inf * [ ...
			sin(theta-aoa_rad), ...
			(cos(theta(1)-aoa_rad) + cos(theta(end)-aoa_rad))]';
	end

	function [cp, cd, cl, cl_KJ] = derive_results( ...
			q, gamma, x_bar, dx, dy, l, theta, D, F, C, G, M_11)
		% DERIVE_RESULTS  Derive relevant results from the H&S solution.
		%
		% Parameters:
		%	q, gamma:
		%	x_bar:
		%	dx, dy:
		%	l, theta:
		%	D, F, C, G, M_11: double(npanel+1, npanel+1)
		%		Internal matrices of the H&S linear system. See
		%		build_system() and theoretical lesson 5 for
		%		implementation details.
		%
		% Returns:
		%	cp:
		%	cd, cl:
		%	cl_KJ:

		% TODO: Try to vectorize v_tan.
		% Tangential velocity on each panel.
		v_tan = zeros(npanel,1);
		% Second value of the tangential velocity.
		val_2 = (D.*F)/2 + C.*G;
		term2 = zeros(npanel, 1);
		for i = 1 : npanel
			val_2(i, i) = 0;
			for j = 1 : npanel
				term2(i) = term2(i) + q(j) * val_2(i, j) / (2*pi);
			end
		end
		% third value of the tangential velocity.
		val_3 = sum(M_11, 2);
		for i = 1 : npanel
			v_tan(i, 1) = ...
				  v_inf * cos(theta(i)-aoa_rad) ...
				- term2(i) ...
				+ gamma/(2*pi) * val_3(i);
		end

		% Pressure coefficient for each panel.
		cp = 1 - (v_tan/v_inf).^2;
		diff_cp = cp(1:end/2) - flip(cp(end/2+1:end));
		diff_pressure = diff_cp * (0.5*rho*v_inf^2);

		% Axial and normal forces.
		Nprime = trapz(x_bar(end/2+1:end), flip(diff_pressure));
		Aprime = trapz(x_bar(end/2+1:end),flip(-diff_pressure'.*(dy(1:end/2)./dx(1:end/2))));

		% Lift and drag forces.
		Lprime = Nprime*cosd(aoa) - Aprime*sind(aoa);
		Dprime = Nprime*sind(aoa) + Aprime*cosd(aoa);

		% Lift and drag coefficients.
		cl = Lprime / (0.5*rho*v_inf^2*c);
		cd = Dprime / (0.5*rho*v_inf^2*c);

		% Lift coefficient from Kutta-Joukowski theorem.
		% To be compared with previously computed cl.
		cl_KJ = 2*gamma*sum(l)/(v_inf*c);
	end

	function plot_cp(x_bar, cp)
		% PLOT_CP  Plot the pressure coefficient distribution along the chord.

		figure('WindowStyle', 'docked');
		hold on;
		
		% Lower surface.
		plot( ...
			x_bar(1:npanel/2)./c, ...
			cp(1:npanel/2), ...
			'color', 'red', ...
			'Linewidth', 1.5);
		% Upper surface.
		plot( ...
			x_bar(npanel/2:npanel)./c, ...
			cp(npanel/2:end), ...
			'color', 'blue', ...
			'Linewidth', 1.5);

		% Dress the plot.
		title('Distribution of the pressure coefficient Cp along the chord');
		xlabel('x/c'); ylabel('Cp');
		set(gca, 'YDir', 'reverse')
		grid;
		legend('Lower surface', 'Upper surface');
	end

	function write_results(X, Y, x_bar, cp, c)
		% WRITE_RESULTS	 Write desired results in external files.

		% Directory holding the result files.
		dir = 'Results/';

		% 1. Airfoil and camber coordinates.

		% Specify the record file name.
		filename = strcat(dir, 'hs-airfoil.csv');

		% Gather the data to store.
		ext_airfoil = [X; Y]';

		% Write in external file.
		writematrix(ext_airfoil, filename);

		% 2. Coefficient of pressure along the airfoil.

		% Specify the record file name.
		filename = strcat( ...
			dir, ...
			'hs-cp', ...
			'-a', num2str(floor(aoa)), ...
			'-v', num2str(floor(v_inf)));

		% Gather the data to store.
		ext_cp_low = [ ...
			x_bar(1:npanel/2)'./c, ...
			cp(1:npanel/2) ];
		ext_cp_up = [ ...
			x_bar(npanel/2:npanel)'./c, ...
			cp(npanel/2:end) ];

		% Write in external file.
		writematrix(ext_cp_up, strcat(filename, '-up.csv'));
		writematrix(ext_cp_low, strcat(filename, '-low.csv'));
	end
end