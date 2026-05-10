%% HELIOCENTRIC TRANSFER OPTIMIZATION: EARTH TO ASTEROID 3908 NYX
% This script computes the optimal transfer trajectory between Earth and 
% Asteroid 3908 Nyx using a hybrid approach: Grid Search for global 
% initialization followed by fmincon (SQP) for local refinement.

clear; clc; close all;

%% 1. PHYSICAL CONSTANTS & CONVERSION FACTORS
AU    = 149597870.7;                % Astronomical Unit [km]
mu    = 1.32712440041279419e+11;    % Solar Gravitational Parameter [km^3/s^2]
R_sun = 696340;                     % Sun Radius [km]
day2sec = 86400;                    % Conversion factor: days to seconds

%% 2. TARGET DATA: ASTEROID 3908 NYX (JPL DATABASE)
a_f  = 1.927925 * AU;               % Semi-major axis [km]
e_f  = 0.459072;                    % Eccentricity [-]
i_f  = deg2rad(2.19);               % Inclination [rad]
OM_f = deg2rad(261.19);             % Longitude of Ascending Node [rad]
om_f = deg2rad(126.66);             % Argument of Perihelion [rad]

%% 3. DEPARTURE DATA: EARTH (INITIAL ORBIT)
a_i  = 1.4946e+08;                  % Semi-major axis [km]
e_i  = 0.016;                       % Eccentricity [-]
i_i  = 9.1920e-05;                  % Inclination [rad]
OM_i = 2.7847;                      % Longitude of Ascending Node [rad]
om_i = 5.2643;                      % Argument of Perihelion [rad]

%% 4. GLOBAL OPTIMIZATION: GRID SEARCH
grid_size = 100;
fprintf('>> Running Grid Search (Size: %d)... ', grid_size);

[th_e_grid, th_n_grid, om_t_grid, dv_grid] = grid_search_s2(...
    a_i, e_i, i_i, OM_i, om_i, ...
    a_f, e_f, i_f, OM_f, om_f, ...
    mu, R_sun, grid_size);

fprintf('DONE.\n');

%% 5. LOCAL OPTIMIZATION: FMINCON (SQP ALGORITHM)
x0 = [th_e_grid, th_n_grid, om_t_grid];  % Initial guess from grid search
lb = [0, 0, 0];                          % Lower bounds for [th_earth, th_nyx, om_tn]
ub = [2*pi, 2*pi, 2*pi];                 % Upper bounds for [th_earth, th_nyx, om_tn]

options = optimoptions('fmincon', ...    % Optimization options
    'Display', 'iter', ...
    'Algorithm', 'sqp', ... 
    'TolFun', 1e-10, ...
    'MaxIterations', 500);

fprintf('\n>> Refining solution with fmincon...\n');
[x_opt, dv_opt] = fmincon(@(x) mission_transfer_nyx(x, a_i, e_i, i_i, ...
    OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu, R_sun), ...
    x0, [], [], [], [], lb, ub, [], options);

%% 6. POST-PROCESSING: RECONSTRUCT OPTIMAL ORBIT
th_earth_opt = x_opt(1);            % Optimal true anomaly for Earth
th_nyx_opt = x_opt(2);              % Optimal true anomaly for Asteroid 3908 Nyx
om_tn_opt = x_opt(3);               % Optimal argument of time for the transfer

[r1_opt, v_earth_opt] = par2car(a_i, e_i, i_i, OM_i, om_i, th_earth_opt, mu);
[r2_opt, v_nyx_opt]   = par2car(a_f, e_f, i_f, OM_f, om_f, th_nyx_opt, mu);

% Compute transfer orbit elements for reporting
h_t_opt  = cross(r1_opt, r2_opt) / norm(cross(r1_opt, r2_opt));         % Specific angular momentum unit vector
i_t_opt  = acos(h_t_opt(3));                                            % Inclination of transfer orbit 
N_t_vect = cross([0; 0; 1], h_t_opt) / norm(cross([0; 0; 1], h_t_opt)); % Node vector for RAAN calculation
OM_t_opt = atan2(N_t_vect(2), N_t_vect(1));                             % RAAN of transfer orbit

if OM_t_opt < 0, OM_t_opt = OM_t_opt + 2*pi; end

% Rotation matrix from Ecliptic to Transfer Plane
RT_ELIO_to_Plane = [1 0 0; 0 cos(i_t_opt) sin(i_t_opt); 0 -sin(i_t_opt) cos(i_t_opt)] * ...
                   [cos(OM_t_opt) sin(OM_t_opt) 0; -sin(OM_t_opt) cos(OM_t_opt) 0; 0 0 1];

r1_p = RT_ELIO_to_Plane * r1_opt;
r2_p = RT_ELIO_to_Plane * r2_opt;

alpha1 = atan2(r1_p(2), r1_p(1));           % True anomaly of r1 in the transfer plane
alpha2 = atan2(r2_p(2), r2_p(1));           % True anomaly of r2 in the transfer plane

while alpha2 < alpha1, alpha2 = alpha2 + 2*pi; end

th_i_t = alpha1 - om_tn_opt;                % True anomaly at departure in the transfer orbit
th_f_t = alpha2 - om_tn_opt;                % True anomaly at arrival in the transfer orbit

r1_n = norm(r1_opt);
r2_n = norm(r2_opt);
e_t  = (r2_n - r1_n) / (r1_n*cos(th_i_t) - r2_n*cos(th_f_t));   % Eccentricity of the transfer orbit
p_t  = r1_n * (1 + e_t*cos(th_i_t));                            % Semi-latus rectum of the transfer orbit
a_t  = p_t / (1 - e_t^2);                                       % Semi-major axis of the transfer orbit    

% Velocity Vectors for DeltaV breakdown
v_coeff = sqrt(mu/p_t);
v1t_pf  = v_coeff * [-sin(th_i_t); e_t + cos(th_i_t); 0];
v2t_pf  = v_coeff * [-sin(th_f_t); e_t + cos(th_f_t); 0];

R_z_om = [cos(om_tn_opt) sin(om_tn_opt) 0; -sin(om_tn_opt) cos(om_tn_opt) 0; 0 0 1];    % Rotation from perifocal to transfer plane
v1t = RT_ELIO_to_Plane' * (R_z_om' * v1t_pf);                                           % Rotate back to Ecliptic frame
v2t = RT_ELIO_to_Plane' * (R_z_om' * v2t_pf);                                           % Rotate back to Ecliptic frame

%% 7. FINAL RESULTS DISPLAY
% Calculate Time of Flight using your TOF_M function
TOF_sec = TOF_M(a_t, e_t, th_i_t, th_f_t, mu);                  
TOF_days = TOF_sec / day2sec;

separator = repmat('=', 1, 55);
fprintf('\n%s\n', separator);
fprintf('          HELIOCENTRIC MISSION DESIGN REPORT\n');
fprintf('%s\n', separator);

fprintf('>> DELTA-V PERFORMANCE:\n');
fprintf('   Grid Search Result:  %10.6f km/s\n', dv_grid);
fprintf('   Fmincon Result:      %10.6f km/s\n', dv_opt);
fprintf('   Optimized Gain:      %10.6f km/s\n', dv_grid - dv_opt);

fprintf('\n>> OPTIMAL TRANSFER ELEMENTS:\n');
fprintf('   Semi-major Axis (a): %12.2f km\n', a_t);
fprintf('   Eccentricity (e):    %12.6f\n', e_t);
fprintf('   Inclination (i):     %12.4f deg\n', rad2deg(i_t_opt));
fprintf('   RAAN (Omega):        %12.4f deg\n', rad2deg(OM_t_opt));
fprintf('   Arg. Pericenter (om):%12.4f deg\n', rad2deg(om_tn_opt));

fprintf('\n>> OPTIMAL PLANET AND ASTEROID ANOMALIES:\n');
fprintf('   True Anomaly at Departure: %12.4f deg\n', rad2deg(th_earth_opt));
fprintf('   True Anomaly at Arrival:   %12.4f deg\n', rad2deg(th_nyx_opt));

fprintf('\n>> MISSION TIMELINE:\n');
fprintf('   Time of Flight:      %12.2f days\n', TOF_days);
fprintf('   Departure Anomaly:   %12.4f deg\n', rad2deg(th_i_t));
fprintf('   Arrival Anomaly:     %12.4f deg\n', rad2deg(th_f_t));

fprintf('\n>> VELOCITY BREAKDOWN (km/s):\n');
fprintf('   V_Earth (Depart):    [%8.4f, %8.4f, %8.4f]\n', v_earth_opt);
fprintf('   V_Trans (Depart):    [%8.4f, %8.4f, %8.4f]\n', v1t);
fprintf('   Delta-V Partenza:     %8.4f km/s\n', norm(v1t - v_earth_opt));
fprintf('   -----------------------------------------------------\n');
fprintf('   V_Nyx (Arrival):     [%8.4f, %8.4f, %8.4f]\n', v_nyx_opt);
fprintf('   V_Trans (Arrival):   [%8.4f, %8.4f, %8.4f]\n', v2t);
fprintf('   Delta-V Arrivo:       %8.4f km/s\n', norm(v_nyx_opt - v2t));
fprintf('%s\n', separator);


