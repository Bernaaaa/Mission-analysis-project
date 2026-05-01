clear 
clc
close all

% N.B. Comment all print statments inside sub functions to avoid cluttering the output during optimization runs !!

% PHISICAL CONSTANTS
mu = 398600.4418; % Earth's gravitational parameter in km^3/s^2
R_earth = 6371; % Earth's radius in km
Lunar_SOI = 384400; % Distance to the Moon in km


% INITIAL ORBITAL PARAMETERS
a = 24400.00; % semi-major axis in km
e = 0.728300; % eccentricity
i = 0.104700; % inclination in radians
OM = 2.361000; % longitude of ascending node in radians
om = 3.107000; % argument of perigee in radians
th = 2.135000; % true anomaly in radians

% FINAL ORBITAL PARAMETERS (derived from final position and velocity)
r_xf = -7090.590200; % final position in km
r_yf = -5612.557300; % final position in km
r_zf = 3948.902900; % final position in km
v_xf = 5.698000; % final velocity in km/s
v_yf = -5.995000; % final velocity in km/s
v_zf = 1.710000; % final velocity in km/s

% Convert final position and velocity to orbital elements
rf = [r_xf; r_yf; r_zf]; % final position vector in km
vf = [v_xf; v_yf; v_zf]; % final velocity vector in km/s

% Convert final position and velocity to orbital elements of the final orbit
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rf, vf, mu);

% Define the search space for optimization
e_min = 0.01; % Minimum eccentricity
e_max = 0.9;  % Maximum eccentricity

rp_min = R_earth + 400; % Minimum perigee radius in km (Earth's radius + 400 km atmospheric buffer)

ra_min = rp_min * (1 + e_max) / (1 - e_max); % Minimum apogee radius in km
ra_max = Lunar_SOI - 400; % Maximum apogee radius in km (Lunar distance minus a buffer to stay within Earth's sphere of influence)

% Create a grid of eccentricity and apogee radius values
e_vec = linspace(e_min, e_max, 100);
ra_vec = linspace(ra_min, ra_max, 100);

% Create a meshgrid for the optimization variables
[E, RA] = meshgrid(e_vec, ra_vec);

% Initialize matrices to store results
DV = zeros(size(E));
DT = zeros(size(E));

% Loop through each combination of eccentricity and apogee radius to calculate the total delta-V and time of flight for the transfer
for k = 1:numel(E)

    % Extract the test eccentricity and apogee radius for the current iteration
    ra_test = RA(k);
    e_test = E(k);
    a_test = ra_test / (1 + e_test);
    rp_test = a_test * (1 - e_test);

    if rp_test < rp_min || ra_test > ra_max
        DV(k) = NaN; % Mark as invalid if perigee is too low or apogee is too high
        DT(k) = NaN; % Mark as invalid if perigee is too low or apogee is too high
        continue
    end

    %bitangent transfer from initial orbit to transfer orbit at perigee
    [DV1a, DV1b, Delta_t_bt1] = bitangentTransfer(a, e, a_test, e_test, 'pa', mu);

    % Time from initial anomaly to perigee of transfer orbit (since it's a bitangent transfer, we can directly use the time from the bitangent transfer function)
    Delta_T1 = Delta_t_bt1; 

    % Time from initial anomaly to perigee
    Delta_T2 = TOF_M(a, e, th, 0, mu); % Time from initial orbit to apogee of transfer orbit

    % Plane change at apogee of transfer orbit
    [DVP, om_p, theta_plane_blt] = changeOrbitalPlane(a_test, e_test, i, OM, om, i_f, OM_f, mu);

    % Time from apogee of transfer orbit to plane change point )
    Delta_T3 = TOF_M(a_test, e_test, pi, theta_plane_blt, mu);

    % Change of argument of pericenter at apogee of transfer orbit
    [DVom, thi_omega_blt, thf_omega_blt] = changePericenterArg(a_test, e_test, om_p, om_f, mu);

    % time from plane change point to pericenter argument change point
    Delta_T4 = TOF_M(a_test, e_test, theta_plane_blt, thi_omega_blt(1), mu);

    % time from pericenter argument change point to apogee (start of bitangent transfer to final orbit)
    Delta_T5 = TOF_M(a_test, e_test, thf_omega_blt(1), pi, mu); 

    %bitangent transfer from transfer orbit to final orbit at apogee
    [DV2a, DV2b, Delta_t_bt2] = bitangentTransfer(a_test, e_test, a_f, e_f, 'ap', mu);

    %time from apogee of transfer orbit to perigee of final orbit (since it's a bitangent transfer, we can directly use the time from the bitangent transfer function)
    Delta_T6 = Delta_t_bt2;

    %time from perigee of final orbit to final anomaly
    Delta_T7 = TOF_M(a_f, e_f, 0, th_f, mu);  
    
    % Total delta-V for the transfer
    DV_tot = abs(DV1a) + abs(DV1b) + abs(DVP) + abs(DVom) + abs(DV2a) + abs(DV2b);
    DV(k) = DV_tot;

    % Total time of flight for the transfer
    DT_tot = Delta_T1 + Delta_T2 + Delta_T3 + Delta_T4 + Delta_T5 + Delta_T6 + Delta_T7;
    DT(k) = DT_tot;
end

%% VISUALIZATION OF THE DATA SET

% Visualize the results for Delta V as function of eccentricity and apogee radius
plotter(E, RA, DV, 'Delta V as function of eccentricity and apogee radius', 'Eccentricity', 'Apogee Radius (km)', 'Delta V (km/s)')

%Visualize the results for Time of Flight as function of eccentricity and apogee radius
plotter(E, RA, DT, 'Time of Flight as function of eccentricity and apogee radius', 'Eccentricity', 'Apogee Radius (km)', 'Time of Flight (s)')

% Visualize the results for Time of Flight as function of eccentricity and apogee radius with logarithmic scale
plotter(E, RA, log10(DT), 'Logarithmic Time of Flight as function of eccentricity and apogee radius', 'Eccentricity', 'Apogee Radius (km)', 'Log10(Time of Flight)')

%% ARRANGE DATA FOR MODEL FITTING

% reshape into arrays
E_flat = E(:);
RA_flat = RA(:);
DT_flat =DT(:); 
DV_flat = DV(:);

% check for Nan values in the data (if any) and handle them before fitting the models
if sum(isnan(DT_flat)) + sum(isnan(DV_flat)) > 0
    valid = ~isnan(DT_flat) & ~isnan(DV_flat);
    E_flat = E_flat(valid);
    RA_flat = RA_flat(valid);
    DT_flat = DT_flat(valid);
    DV_flat = DV_flat(valid);
end

%% MODEL FITTING

%POLYNOMIAL FITTING
% Fit a polynomial surface to the valid data points
ft_DT= fit([E_flat, RA_flat], DT_flat, 'poly12'); % 'poly12' means a polynomial of degree 1 in E and degree 2 in RA (total degree 3)
ft_DV = fit([E_flat, RA_flat], DV_flat, 'poly23'); % 'poly23' means a polynomial of degree 2 in E and degree 3 in RA (total degree 5)

% VISUALIZATION OF POLYNOMIAL FITTING
% Reshape t_val back onto the full 2D grid
DT_poly_grid = reshape(ft_DT(E(:), RA(:)), size(E));
DV_poly_grid = reshape(ft_DV(E(:), RA(:)), size(E));

% Visualize the fitted polynomial surface for Time of Flight
plotter(E, RA, DT_poly_grid, 'Fitted Polynomial(poly12) Time of Flight Surface', 'Eccentricity', 'Apogee Radius (km)', 'Time of Flight (s)')
plotter(E, RA, DV_poly_grid, 'Fitted Polynomial(poly23) Delta V Surface', 'Eccentricity', 'Apogee Radius (km)', 'Delta V (km/s)')

% SPLINE FITTING
% Fit a scatteredInterpolant spline
F_DT = scatteredInterpolant(E_flat, RA_flat, DT_flat, 'natural', 'none');
F_DV = scatteredInterpolant(E_flat, RA_flat, DV_flat, 'natural', 'none');

% VISUALIZATION OF SCATTERED INTERPOLANT SPLINE FITTING
% Evaluate the spline at the grid points
DT_spline2 = F_DT(E, RA);
DV_spline2 = F_DV(E, RA);

% Visualize the scatteredInterpolant surface for Time of Flight
plotter(E, RA, DT_spline2, 'ScatteredInterpolant (natural) Time of Flight Surface', 'Eccentricity', 'Apogee Radius (km)', 'Time of Flight (s)')
plotter(E, RA, DV_spline2, 'ScatteredInterpolant (natural) Delta V Surface', 'Eccentricity', 'Apogee Radius (km)', 'Delta V (km/s)')

% CUBIC SPLINE FITTING
% cubic interpolation (this will also return NaN for points outside the convex hull of valid data)
DT_spline = griddata(E_flat, RA_flat, DT_flat, E, RA, 'cubic');
DV_spline = griddata(E_flat, RA_flat, DV_flat, E, RA, 'cubic');

%VISUALIZATION OF CUBIC SPLINE FITTING
% Visualize the cubic interpolated surface for Time of Flight
plotter(E, RA, DT_spline, 'Cubic Interpolated Time of Flight Surface', 'Eccentricity', 'Apogee Radius (km)', 'Time of Flight (s)')
plotter(E, RA, DV_spline, 'Cubic Interpolated Delta V Surface', 'Eccentricity', 'Apogee Radius (km)', 'Delta V (km/s)')
% calculate metrics for each model (RMSE and R²)

%% MODEL EVALUATION

% Metrics for Time of Flight
[rmse_DT_poly,   r2_DT_poly]   = model_metrics(DT_poly_grid, DT);
[rmse_DT_spline, r2_DT_spline] = model_metrics(DT_spline,    DT);
[rmse_DT_spline2, r2_DT_spline2] = model_metrics(DT_spline2, DT);

% Metrics for Delta V
[rmse_DV_poly,   r2_DV_poly]   = model_metrics(DV_poly_grid, DV);
[rmse_DV_spline, r2_DV_spline] = model_metrics(DV_spline,    DV);
[rmse_DV_spline2, r2_DV_spline2] = model_metrics(DV_spline2, DV);


% Display results in a formatted table
fprintf('\n%-20s %10s %10s\n', 'Modello', 'RMSE', 'R²');
fprintf('%s\n', repmat('-', 1, 42));

fprintf('%-20s %10.4f %10.6f\n', 'DT  poly31',  rmse_DT_poly,    r2_DT_poly);
fprintf('%-20s %10.4f %10.6f\n', 'DT  spline',  rmse_DT_spline,  r2_DT_spline);
fprintf('%-20s %10.4f %10.6f\n', 'DT  spline2', rmse_DT_spline2, r2_DT_spline2);

fprintf('%-20s %10.4f %10.6f\n', 'DV  poly31',  rmse_DV_poly,    r2_DV_poly);
fprintf('%-20s %10.4f %10.6f\n', 'DV  spline',  rmse_DV_spline,  r2_DV_spline);
fprintf('%-20s %10.4f %10.6f\n', 'DV  spline2', rmse_DV_spline2, r2_DV_spline2);

%% OPTIMIZATION


% normlizated weighted optimization
% obj = (1 - w) * F_DT + w * F_DV
% w = 0   → minimize DT only
% w = 1   → minimize DV only
% w = 0.5 → equal weight to both DT and DV

w = 0.5;

% normalization factors (use mean or max of valid data to scale the objectives) 
DT_ref = mean(DT_flat); 
DV_ref = mean(DV_flat);

%OBJECTIVE FUNCTION FOR OPTIMIZATION
% combined objective function (normalized and weighted)
obj = @(e, ra) (1 - w) * F_DT(e, ra) ./ DT_ref + ...
        w  * F_DV(e, ra) ./ DV_ref;

% Define the objective function for optimization (wrap the combined objective in a function that takes a vector input)
obj_vec = @(x) obj(x(1), x(2));

% Set bounds for optimization
lb = [min(E_flat),  min(RA_flat)];
ub = [max(E_flat),  max(RA_flat)];

% Optimization options
options = optimoptions('fmincon', ...
    'Display',       'off', ...
    'Algorithm',     'sqp', ...
    'TolFun',        1e-8, ...
    'TolX',          1e-8, ...
    'MaxIterations', 500);

%multiple local minima analysis (optional, requires global optimization toolbox)
n_start = 50; % numero di punti di partenza casuali
results = zeros(n_start, 3); % [e_opt, ra_opt, f_opt]

% Genera punti di partenza casuali dentro il dominio valido
rng(42) % riproducibilità

e_starts  = e_min + (e_max  - e_min)  * rand(n_start, 1);
ra_starts = ra_min + (ra_max - ra_min) * rand(n_start, 1);

for k = 1:n_start
    x0_k = [e_starts(k), ra_starts(k)];
    
    % Salta se il punto di partenza è fuori dal convex hull (NaN)
    if isnan(F_DT(x0_k(1), x0_k(2))) || isnan(F_DV(x0_k(1), x0_k(2)))
        results(k, :) = [NaN, NaN, Inf];
        continue
    end
    
    try
        [xk, fk] = fmincon(obj_vec, x0_k, [], [], [], [], lb, ub, [], options);
        results(k, :) = [xk(1), xk(2), fk];
    catch
        results(k, :) = [NaN, NaN, Inf];
    end
end

% Prendi il migliore tra tutti i run validi

valid_runs = results(:,3) < Inf & ~isnan(results(:,3));
[ft_opt_ms, idx_best] = min(results(valid_runs, 3));
best = results(valid_runs, :);
best = best(idx_best, :);

e_opt_ms  = best(1);
ra_opt_ms = best(2);
DT_opt_ms = F_DT(e_opt_ms, ra_opt_ms);
DV_opt_ms = F_DV(e_opt_ms, ra_opt_ms);

fprintf('\n--- Multi-start fmincon ---\n')
fprintf('e  : %.4f\n',       e_opt_ms)
fprintf('r_a: %.2f km\n',    ra_opt_ms)
fprintf('DT_opt_ms : %.2f s (%.2f ore)\n', DT_opt_ms, DT_opt_ms/3600)
fprintf('DV_opt_ms : %.4f km/s\n',   DV_opt_ms)

% Plot the combined objective function surface
plotter(E, RA, obj(E, RA), sprintf('Funzione obiettivo combinata (w = %.2f)', w), 'Eccentricity', 'Apogee Radius (km)', 'Funzione obiettivo')
hold on
plot3(e_opt_ms, ra_opt_ms, ft_opt_ms, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r')
xlabel('e'); ylabel('r_a (km)'); zlabel('Funzione obiettivo')
title(sprintf('Funzione obiettivo combinata (w = %.2f) con multi-start', w))
legend('Superficie', 'Ottimo')
hold off

e_test = e_opt_ms;
ra_test = ra_opt_ms;
a_test = ra_test / (1 + e_test);
rp_test = a_test * (1 - e_test);

%bitangent transfer from initial orbit to transfer orbit at perigee
[DV1a, DV1b, Delta_t_bt1] = bitangentTransfer(a, e, a_test, e_test, 'pa', mu);
Delta_T1 = TOF_M(a, e, th, 0, mu); % Time from initial orbit to apogee of transfer orbit

% Plane change at apogee of transfer orbit
[DVP, om_p, theta_plane_blt] = changeOrbitalPlane(a_test, e_test, i, OM, om, i_f, OM_f, mu);
Delta_T2 = TOF_M(a_test, e_test, pi, theta_plane_blt, mu);

    % Change of argument of pericenter at apogee of transfer orbit
[DVom, thi_omega_blt, thf_omega_blt] = changePericenterArg(a_test, e_test, om_p, om_f, mu);
Delta_T3 = TOF_M(a_test, e_test, theta_plane_blt, thi_omega_blt(1), mu); % Time from plane change point to pericenter argument change point
Delta_T4 = TOF_M(a_test, e_test, thf_omega_blt(1), pi, mu); % Time from pericenter argument change point to final orbit

%bitangent transfer from transfer orbit to final orbit at apogee
[DV2a, DV2b, Delta_t_bt2] = bitangentTransfer(a_test, e_test, a_f, e_f, 'ap', mu);
Delta_T5 = TOF_M(a_f, e_f, 0, th_f, mu); % Time from pericenter argument change point to final orbit 
    
% Total delta-V for the transfer
DV_tot = abs(DV1a) + abs(DV1b) + abs(DVP) + abs(DVom) + abs(DV2a) + abs(DV2b);


% Total time of flight for the transfer
DT_tot = Delta_T1 + Delta_T2 + Delta_T3 + Delta_T4 + Delta_T5 + Delta_t_bt1 + Delta_t_bt2;

fprintf('\n--- Soluzione ottima ESATTA ---\n')
fprintf('Time of Flight totale: %.2f s (%.2f ore)\n', DT_tot, DT_tot/3600)
fprintf('Delta V totale: %.4f km/s\n', DV_tot)



options = optimoptions('patternsearch', ...
    'Display','off', ...
    'MaxIterations',300);

for k = 1:n_start
    x0_k = [e_starts(k), ra_starts(k)];
    [x_opt, f_opt] = patternsearch(obj_vec, x0_k, [], [], [], [], lb, ub, options);
    results(k, :) = [x_opt(1), x_opt(2), f_opt];
end

valid_runs = results(:,3) < Inf & ~isnan(results(:,3));
[ft_opt_ms, idx_best] = min(results(valid_runs, 3));
best = results(valid_runs, :);
best = best(idx_best, :);

e_opt_ms  = best(1);
ra_opt_ms = best(2);
DT_opt_ms = F_DT(e_opt_ms, ra_opt_ms);
DV_opt_ms = F_DV(e_opt_ms, ra_opt_ms);

fprintf('\n--- Multi-start patternsearch ---\n')
fprintf('e  : %.4f\n',       e_opt_ms)
fprintf('r_a: %.2f km\n',    ra_opt_ms)
fprintf('DT_opt_ms : %.2f s (%.2f ore)\n', DT_opt_ms, DT_opt_ms/3600)
fprintf('DV_opt_ms : %.4f km/s\n',   DV_opt_ms)

