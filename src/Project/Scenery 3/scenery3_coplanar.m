% =========================================================================
% INTERPLANETARY MISSION ANALYSIS - SCENARIO 3: PATCHED CONICS
% Calculation of Departure (Earth) and Arrival (Nyx) Hyperbolic Trajectories
% =========================================================================

clear; clc; close all;

%% 0. PATH SETUP
% Ensure the common functions folder is included in the MATLAB path
currentDir = fileparts(mfilename('fullpath'));
funcFolder = fullfile(currentDir, '..', '..', 'lab');

if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
else
    warning('Warning: Functions folder not found at: %s', funcFolder);
end

%% 1. PHYSICAL CONSTANTS & CONVERSION FACTORS
% -------------------------------------------------------------------------
% Fundamental masses and dimensions
M_sun   = 1.989e30;       % [kg] Mass of the Sun
M_earth = 5.972e24;       % [kg] Mass of the Earth
M_nyx   = 1.0472e12;      % [kg] Mass of Asteroid 3908 Nyx
D_nyx   = 1.0;            % [km] Diameter of 3908 Nyx
AU      = 149597870.7;    % [km] Astronomical Unit
R_earth = 6378.137;       % [km] Earth's equatorial radius

% Gravitational parameters (Converted to km-based units: km^3 / (kg * s^2))
G        = 6.67430e-20;   
mu_earth = G * M_earth;   % [km^3/s^2] Earth's gravitational parameter
mu_sun   = G * M_sun;     % [km^3/s^2] Sun's gravitational parameter
mu_nyx   = G * M_nyx;     % [km^3/s^2] Nyx's gravitational parameter


%% 2. ORBITAL PARAMETERS & POSITIONS
% -------------------------------------------------------------------------
% --- Earth Parking Orbit State Vector (From Scenery 1) ---
% Position [km] and Velocity [km/s] relative to Earth
ri = [-7090.590200; -5612.557300; 3948.902900];
vi = [5.698000; -5.995000; 1.710000];

% Convert Cartesian state vector to Keplerian elements
[a_pi, e_pi, i_pi, OM_pi, om_pi, th_pi] = car2par(ri, vi, mu_earth);

% --- Heliocentric Orbits (From Scenery 2) ---
% Earth's heliocentric parameters
a_earth  = 1.4946e+08;           % [km] Semi-major axis
e_earth  = 0.016;                % [-] Eccentricity
th_earth = deg2rad(252.1890);    % [rad] True Anomaly at departure

% Nyx's heliocentric parameters
a_nyx  = 1.927925 * AU;          % [km] Semi-major axis
e_nyx  = 0.459072;               % [-] Eccentricity
th_nyx = deg2rad(28.4244);       % [rad] True Anomaly at arrival


%% 3. SPHERES OF INFLUENCE (SOI) CALCULATION
% -------------------------------------------------------------------------
% Heliocentric distance of Earth at departure
p_earth = a_earth * (1 - e_earth^2);   
R_earth_departure = p_earth / (1 + e_earth * cos(th_earth));

% Earth's Sphere of Influence
earth_SOI = R_earth_departure * (M_earth / M_sun)^(2/5);  

% Heliocentric distance of Nyx at arrival
p_nyx = a_nyx * (1 - e_nyx^2);         
R_nyx_arrival = p_nyx / (1 + e_nyx * cos(th_nyx));

% Nyx's Sphere of Influence
nyx_SOI = R_nyx_arrival * (M_nyx / M_sun)^(2/5);  


%% 4. ARRIVAL HYPERBOLA (ASTEROID NYX CAPTURE)
% -------------------------------------------------------------------------
% Based on Scenery 2, the hyperbolic excess velocity at arrival is known
v_inf_nyx = 4.0173; % [km/s]

% Hyperbola semi-major axis
a_h_nyx = -mu_nyx / (v_inf_nyx^2);

% Define search grid for periapsis radius (From surface to SOI edge)
r_ph_grid = linspace((D_nyx / 2) + .1, nyx_SOI, 1000);
best_dv_nyx = inf;

% Find the optimal periapsis radius that minimizes capture Delta-V
for rp = r_ph_grid
    
    % Circular orbit velocity at current altitude
    v_ocn = sqrt(mu_nyx / rp);
    
    % Hyperbola geometry
    eh = 1 - (rp / a_h_nyx);
    
    % Velocity at hyperbolic periapsis 
    vph = sqrt(v_inf_nyx^2 + 2 * mu_nyx / rp);
    
    % Required Delta-V for capture maneuver
    deltaV = abs(v_ocn - vph);
    
    if deltaV < best_dv_nyx
        best_dv_nyx = deltaV;
        best_rp_nyx = rp;
        best_eh_nyx = eh;
        best_impact_n_nyx = -a_h_nyx * sqrt(eh^2 - 1);
        best_deflection_nyx = 2 * asin(1 / eh);
    end
end 


%% 5. DEPARTURE HYPERBOLA (EARTH ESCAPE)
% -------------------------------------------------------------------------
% Based on Scenery 2, the hyperbolic excess velocity at departure is known
v_inf_earth = 2.9986; % [km/s]

% Hyperbola semi-major axis
a_h_earth = -mu_earth / (v_inf_earth^2);

% The optimal maneuver occurs at the periapsis of the parking orbit (Oberth effect)
r_pi_earth = a_pi * (1 - e_pi); 
v_pi_earth = sqrt(mu_earth / r_pi_earth) * sqrt(1 + e_pi); % Real velocity at periapsis

% Hyperbola geometry
e_h_earth = 1 - (r_pi_earth / a_h_earth); 

% Velocity required at hyperbolic periapsis to escape with v_inf_earth
vph_earth = sqrt(v_inf_earth^2 + 2 * mu_earth / r_pi_earth);

% Required Delta-V for departure maneuver
deltaV_earth = abs(vph_earth - v_pi_earth);

% Calculate geometric parameters
impact_n_earth = -a_h_earth * sqrt(e_h_earth^2 - 1);
deflection_earth = 2 * asin(1 / e_h_earth);


%% 6. MISSION REPORT OUTPUT
% -------------------------------------------------------------------------
fprintf('\n=========================================================\n');
fprintf('           SCENARIO 3: PATCHED CONICS ANALYSIS           \n');
fprintf('=========================================================\n\n');

fprintf('--- SPHERES OF INFLUENCE (SOI) ---\n');
fprintf('Earth SOI at departure : %10.2f km\n', earth_SOI);
fprintf('Nyx SOI at arrival     : %10.4f km\n\n', nyx_SOI);

fprintf('--- PHASE 1: EARTH DEPARTURE ---\n');
fprintf('Maneuver Delta-V       : %10.4f km/s\n', deltaV_earth);
fprintf('Maneuver Altitude (rp) : %10.4f km\n', r_pi_earth);
fprintf('Hyperbola Eccentricity : %10.4f\n', e_h_earth);
fprintf('Deflection Angle       : %10.4f deg\n', rad2deg(deflection_earth));
fprintf('Impact Parameter       : %10.4f km\n\n', impact_n_earth);

fprintf('--- PHASE 3: NYX ARRIVAL (CAPTURE) ---\n');
fprintf('Minimum Delta-V        : %10.4f km/s\n', best_dv_nyx);
fprintf('Optimal Altitude (rp)  : %10.4f km\n', best_rp_nyx);
fprintf('Hyperbola Eccentricity : %10.4e\n', best_eh_nyx);
fprintf('Deflection Angle       : %10.4e deg\n', rad2deg(best_deflection_nyx));
fprintf('Impact Parameter       : %10.4f km\n', best_impact_n_nyx);

fprintf('\n=========================================================\n');
fprintf('TOTAL PATCHED CONICS DELTA-V: %.4f km/s\n', abs(deltaV_earth) + abs(best_dv_nyx));
fprintf('=========================================================\n');

%% 7. 3D VISUALIZATION
% -------------------------------------------------------------------------

earth_plot_data = struct(...
    'mu', mu_earth, ...
    'R', R_earth, ...
    'a_park', a_pi, ...
    'e_park', e_pi, ....
    'i', i_pi, ...
    'OM', OM_pi, ...
    'om', om_pi, ...
    'a_hyp', a_h_earth, ...
    'e_hyp', e_h_earth ...
);

nyx_plot_data = struct(...
    'mu', mu_nyx, ...
    'R', D_nyx/2, ...
    'r_circ', best_rp_nyx, ...
    'a_hyp', a_h_nyx, ...
    'e_hyp', best_eh_nyx ...
);

scenery3_plot_coplanar(earth_plot_data, nyx_plot_data);