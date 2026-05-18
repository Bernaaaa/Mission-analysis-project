% =========================================================================
% INTERPLANETARY MISSION ANALYSIS - SCENARIO 3: PATCHED CONICS
% Calculation of Departure (Earth) and Arrival (Nyx) Hyperbolic Trajectories
% =========================================================================

clear; clc; close all;

%% 0. PATH SETUP
% -------------------------------------------------------------------------
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

% Gravitational parameters [km^3/s^2]
G        = 6.67430e-20;   
mu_earth = G * M_earth;   
mu_sun   = G * M_sun;     
mu_nyx   = G * M_nyx;     


%% 2. ORBITAL PARAMETERS & VECTORS
% -------------------------------------------------------------------------
% --- Earth Parking Orbit (Scenario 1 - ECI Frame) ---
ri = [-7090.590200; -5612.557300; 3948.902900];
vi = [5.698000; -5.995000; 1.710000];
[a_pi, e_pi, i_pi, OM_pi, om_pi, th_pi] = car2par(ri, vi, mu_earth);

% --- Heliocentric Orbits (Scenario 2 - Ecliptic Frame) ---
a_earth  = 1.4946e+08;           % [km] Semi-major axis Earth
e_earth  = 0.016;                % [-] Eccentricity Earth
th_earth = deg2rad(252.1890);    % [rad] True Anomaly Earth at departure

a_nyx  = 1.927925 * AU;          % [km] Semi-major axis Nyx
e_nyx  = 0.459072;               % [-] Eccentricity Nyx
th_nyx = deg2rad(28.4244);       % [rad] True Anomaly Nyx at arrival

% --- Coordinate Transformation (Ecliptic to ECI) ---
epsilon = deg2rad(23.45);        % [rad] Obliquity of the ecliptic
T_ecliptic_to_ECI = [1, 0, 0; 
                     0, cos(epsilon), -sin(epsilon); 
                     0, sin(epsilon),  cos(epsilon)];

% --- Earth Departure (V_inf calculation) ---
V_Earth        = [2.9768; 29.5101; -0.0026];
V_Trans_Depart = [3.1275; 32.4458;  0.5897];

v_inf_vec_earth = V_Trans_Depart - V_Earth; 
v_inf_vec_earth_eci = T_ecliptic_to_ECI * v_inf_vec_earth; 
v_inf_earth = norm(v_inf_vec_earth_eci); 

% --- Nyx Arrival (V_inf calculation) ---
V_Nyx           = [-25.2398;  23.2036; -1.0897];
V_Trans_Arrival = [-22.5634;  20.5617;  0.3229];

v_inf_vec_nyx = V_Trans_Arrival - V_Nyx; % Flipped to target reference
v_inf_nyx = norm(v_inf_vec_nyx); 


%% 3. SPHERES OF INFLUENCE (SOI) 
% -------------------------------------------------------------------------
% Earth SOI
p_earth = a_earth * (1 - e_earth^2);   
R_earth_departure = p_earth / (1 + e_earth * cos(th_earth));
earth_SOI = R_earth_departure * (M_earth / M_sun)^(2/5);  

% Nyx SOI
p_nyx = a_nyx * (1 - e_nyx^2);         
R_nyx_arrival = p_nyx / (1 + e_nyx * cos(th_nyx));
nyx_SOI = R_nyx_arrival * (M_nyx / M_sun)^(2/5);  


%% 4. ARRIVAL HYPERBOLA (ASTEROID NYX CAPTURE)
% -------------------------------------------------------------------------
a_h_nyx = -mu_nyx / (v_inf_nyx^2);

r_ph_grid = linspace(D_nyx/2 + 0.1, nyx_SOI, 10000); 
best_dv_nyx = inf;

for rp = r_ph_grid
    v_ocn = sqrt(mu_nyx / rp);
    eh = 1 - (rp / a_h_nyx);
    vph = sqrt(v_inf_nyx^2 + 2 * mu_nyx / rp);
    
    deltaV = abs(v_ocn - vph);
    
    if deltaV < best_dv_nyx
        best_dv_nyx         = deltaV;
        best_rp_nyx         = rp;
        best_eh_nyx         = eh;
        best_impact_n_nyx   = -a_h_nyx * sqrt(eh^2 - 1);
        best_deflection_nyx = 2 * asin(1 / eh);
    end
end 


%% 5. DEPARTURE HYPERBOLA (EARTH ESCAPE 3D)
% -------------------------------------------------------------------------
a_h_earth = -mu_earth / (v_inf_earth^2);

% Orbital properties at periapsis (Oberth effect)
r_pi_earth = a_pi * (1 - e_pi); 
v_pi_earth = sqrt(mu_earth / r_pi_earth) * sqrt(1 + e_pi); % Velocity of parking orbit

% Velocity required at hyperbolic periapsis 
vph_earth = sqrt(v_inf_earth^2 + 2 * mu_earth / r_pi_earth); % Hyperbolic velocity

% --- 3D Vectorial Delta-V Calculation (Carnot Theorem) ---
% 1. Angular momentum vectors (defining the two orbital planes)
h_park = cross(ri, vi);
h_hyperbola = cross(ri, v_inf_vec_earth_eci); 

% 2. Plane change angle (misalignment)
cos_Delta_i = dot(h_park, h_hyperbola) / (norm(h_park) * norm(h_hyperbola));
Delta_i = acos(cos_Delta_i); 

% 3. Total Delta-V (Magnitude of the vectorial difference)
deltaV_earth = sqrt(v_pi_earth^2 + vph_earth^2 - 2 * v_pi_earth * vph_earth * cos(Delta_i));
% ---------------------------------------------------------

% Geometric parameters
e_h_earth        = 1 - (r_pi_earth / a_h_earth); 
impact_n_earth   = -a_h_earth * sqrt(e_h_earth^2 - 1);
deflection_earth = 2 * asin(1 / e_h_earth);


%% 6. MISSION REPORT OUTPUT
% -------------------------------------------------------------------------
fprintf('\n=========================================================\n');
fprintf('       PATCHED CONICS MULTI-FRAME ANALYSIS REPORT        \n');
fprintf('=========================================================\n\n');

fprintf('--- SPHERES OF INFLUENCE ---\n');
fprintf(' Earth SOI Radius : %10.2f km\n', earth_SOI);
fprintf(' Nyx SOI Radius   : %10.2f km\n\n', nyx_SOI);

fprintf('--- DEPARTURE (EARTH ESCAPE 3D) ---\n');
fprintf(' V_inf Magnitude  : %10.4f km/s\n', v_inf_earth);
fprintf(' Plane Sfasamento : %10.2f degrees\n', rad2deg(Delta_i));
fprintf(' Optimal Periapsis: %10.2f km\n', r_pi_earth);
fprintf(' Eccentricity     : %10.4f\n', e_h_earth);
fprintf(' Impact Parameter : %10.2f km\n', impact_n_earth);
fprintf(' Deflection Angle : %10.2f degrees\n', rad2deg(deflection_earth));
fprintf(' --------------------------------------\n');
fprintf(' TOTAL DELTA-V    : %10.4f km/s\n\n', deltaV_earth);

fprintf('--- ARRIVAL (NYX CAPTURE) ---\n');
fprintf(' V_inf Magnitude  : %10.4f km/s\n', v_inf_nyx);
fprintf(' Optimal Periapsis: %10.2f km\n', best_rp_nyx);
fprintf(' Eccentricity     : %10.4e\n', best_eh_nyx);
fprintf(' Impact Parameter : %10.2f km\n', best_impact_n_nyx);
fprintf(' Deflection Angle : %10.2e degrees\n', rad2deg(best_deflection_nyx));
fprintf(' --------------------------------------\n');
fprintf(' TOTAL DELTA-V    : %10.4f km/s\n', best_dv_nyx);
fprintf('\n=========================================================\n');