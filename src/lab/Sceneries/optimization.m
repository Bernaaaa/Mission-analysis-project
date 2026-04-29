clear 
clc
close all

% Define constants and initial conditions
mu = 398600.4418;

% Initial orbital parameters of the orbit
a = 24400.00;
e = 0.728300;
i = 0.104700;
OM = 2.361000;
om = 3.107000;
th = 2.135000;

% Final orbital parameters of the orbit
r_xf = -7090.590200;
r_yf = -5612.557300;
r_zf = 3948.902900;
v_xf = 5.698000;
v_yf = -5.995000;
v_zf = 1.710000;

% Convert final position and velocity to orbital elements
rf = [r_xf; r_yf; r_zf];
vf = [v_xf; v_yf; v_zf];

% Convert initial Cartesian coordinates to orbital elements
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rf, vf, mu);

% Define the search space for optimization
e_min = 0.01; % Minimum eccentricity
e_max = 0.9;  % Maximum eccentricity
rp_min = 6371 + 400; % Minimum perigee radius in km (Earth's radius + 400 km atmospheric buffer)
ra_min = rp_min * (1 + e_max) / (1 - e_max); % Minimum apogee radius in km
ra_max = 350000; % Maximum apogee radius in km

% Create a grid of eccentricity and apogee radius values
e_vec = linspace(e_min, e_max, 100);
ra_vec = linspace(ra_min, ra_max, 100);

[E, RA] = meshgrid(e_vec, ra_vec);

% Initialize matrices to store results
DV = zeros(size(E));
DT = zeros(size(E));



for k = 1:numel(E)

    % Extract the test eccentricity and apogee radius for the current iteration
    ra_test = RA(k);
    e_test = E(k);
    a_test = ra_test / (1 + e_test);
    rp_test = a_test * (1 - e_test);

    %check if the perigee of the transfer orbit is above the Earth's atmosphere (assumed to be 200 km above Earth's surface, so 6578 km from Earth's center)
    if rp_test < 6371 + 400
        DV(k) = NaN; % Mark as invalid
        DT(k) = NaN; % Mark as invalid
        continue;
    end

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
    DV(k) = DV_tot;

    % Total time of flight for the transfer
    DT_tot = Delta_T1 + Delta_T2 + Delta_T3 + Delta_T4 + Delta_T5 + Delta_t_bt1 + Delta_t_bt2;
    DT(k) = DT_tot;
end

% Visualize the results for Delta V as function of eccentricity and apogee radius
plotter(E, RA, DV, 'Delta V as function of eccentricity and apogee radius')

%Visualize the results for Time of Flight as function of eccentricity and apogee radius
plotter(E, RA, DT, 'Time of Flight as function of eccentricity and apogee radius')

% Visualize the results for Time of Flight as function of eccentricity and apogee radius with logarithmic scale
plotter(E, RA, log10(DT), 'Logarithmic Time of Flight as function of eccentricity and apogee radius')

% reshape into arrays
E_flat = E(:);
RA_flat = RA(:);
DT_flat =DT(:); 
DV_flat = DV(:);

% Check if all invalid points are due to perigee being below the atmospheric threshold
if sum(isnan(DT_flat)) + sum(isnan(DV_flat)) > 0
    valid = ~isnan(DT_flat) & ~isnan(DV_flat);
    E_flat = E_flat(valid);
    RA_flat = RA_flat(valid);
    DT_flat = DT_flat(valid);
    DV_flat = DV_flat(valid);
end

% Fit a polynomial surface to the valid data points (e.g., using poly11 for a simple linear fit or poly23 for a more complex fit)
ft_DT= fit([E_flat, RA_flat], DT_flat, 'poly12');
ft_DV = fit([E_flat, RA_flat], DV_flat, 'poly23');

% Reshape t_val back onto the full 2D grid
DT_poly_grid = reshape(ft_DT(E(:), RA(:)), size(E));
DV_poly_grid = reshape(ft_DV(E(:), RA(:)), size(E));

% Visualize the fitted polynomial surface for Time of Flight
plotter(E, RA, DT_poly_grid, 'Fitted Polynomial(poly12) Time of Flight Surface')

% Fit a scatteredInterpolant spline surface to the valid data points (natural neighbor interpolation)
F_DT = scatteredInterpolant(E_flat, RA_flat, DT_flat, 'natural', 'none');
F_DV = scatteredInterpolant(E_flat, RA_flat, DV_flat, 'natural', 'none');

% Evaluate the spline at the grid points
DT_spline2 = F_DT(E, RA);
DV_spline2 = F_DV(E, RA);

% Visualize the scatteredInterpolant surface for Time of Flight
plotter(E, RA, DT_spline2, 'ScatteredInterpolant (natural) Time of Flight Surface')

% cubic interpolation (this will also return NaN for points outside the convex hull of valid data)
DT_spline = griddata(E_flat, RA_flat, DT_flat, E, RA, 'cubic');
DV_spline = griddata(E_flat, RA_flat, DV_flat, E, RA, 'cubic');

% Visualize the cubic interpolated surface for Time of Flight
plotter(E, RA, DT_spline, 'Cubic Interpolated Time of Flight Surface')
% calculate metrics for each model (RMSE and R²)

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


% normlizated weighted optimization
% obj = (1 - w) * F_DT + w * F_DV
% w = 0   → minimize DT only
% w = 1   → minimize DV only
% w = 0.5 → equal weight to both DT and DV

w = 0.5;

% normalization factors (use mean or max of valid data to scale the objectives) 
DT_ref = mean(DT_flat); 
DV_ref = mean(DV_flat);


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
plotter(E, RA, obj(E, RA), sprintf('Funzione obiettivo combinata (w = %.2f)', w))
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

