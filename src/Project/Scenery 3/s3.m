clear; clc; close all;

%% --- SETUP PERCORSI UNIVERSALE ---
currentDir = fileparts(mfilename('fullpath'));

funcFolder = fullfile(currentDir, '..', '..', 'lab');

if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
else
    warning('Attenzione: La cartella delle funzioni non è stata trovata in: %s', funcFolder);
end

%% 1. PHYSICAL CONSTANTS & CONVERSION FACTORS

% --- Physical Constants ---
M_sun = 1.989e30;                       % Mass of the Sun in kg
M_earth = 5.972e24;                     % Mass of the Earth in kg
M_nyx = 1.0472e12;                      % Mass of 3908 Nyx in kg
D_nyx = 1.0;                            % Diameter of 3908 Nyx in km
AU = 149597870.7;                       % Astronomical Unit in km
G = 6.67430e-20;                        % km^3 kg^-1 s^-2                        
mu_earth = G * M_earth;                   % Gravitational parameter of Earth in km^3/s^2
mu_sun = G * M_sun;                       % Gravitational parameter of Sun in km^3/s^2
mu_nyx = G * M_nyx;                       % Gravitational parameter of Nyx in km^3/s^2

%% 2. ORBITAL PARAMETERS & POSITIONS

% --- Earth parking orbit (Scenery 1) ---
r_xpi = -7090.590200;
r_ypi = -5612.557300;
r_zpi = 3948.902900;
v_xpi = 5.698000;
v_ypi = -5.995000;
v_zpi = 1.710000;

ri = [r_xpi; r_ypi; r_zpi];
vi = [v_xpi; v_ypi; v_zpi];

[a_pi, e_pi, i_pi, OM_pi, om_pi, th_pi] = car2par(ri, vi, mu_earth);

% --- Earth orbit (Scenery 2) ---
a_earth  = 1.4946e+08;                  % Semi-major axis [km]
e_earth  = 0.016;                       % Eccentricity [-]
i_earth  = 9.1920e-05;                  % Inclination [rad]
OM_earth = 2.7847;                      % Longitude of Ascending Node [rad]
om_earth = 5.2643;                      % Argument of Perihelion [rad]
th_earth = deg2rad(252.1890);           % True Anomaly [rad]

% --- Nyx orbit (Scenery 2) ---
a_nyx  = 1.927925 * AU;               % Semi-major axis [km]
e_nyx  = 0.459072;                    % Eccentricity [-]
i_nyx  = deg2rad(2.19);               % Inclination [rad]
OM_nyx = deg2rad(261.19);             % Longitude of Ascending Node [rad]
om_nyx = deg2rad(126.66);             % Argument of Perihelion [rad]
th_nyx = deg2rad(28.4244);            % True Anomaly [rad]



%% Strategy
% nyx hyperbolic
a_h = -mu_nyx/(4.0173^2);
r_ph = linspace(D_nyx/2 + 0.01, 2000, 100000);
best_dv = inf;

for rp = r_ph
    v_ocn = sqrt(mu_nyx/rp);
    eh = 1 - rp/a_h;

    vph = sqrt(4.0173^2 + 2*mu_nyx/rp);
    impact_n = -a_h*sqrt(eh^2-1);
    deltaV = abs(v_ocn - vph);

    if deltaV < best_dv
        best_dv = deltaV;
        best_rp = rp;
    end
end 

fprintf('Nyx deltaV: %.4f km/s at rp = %.4f km\n', best_dv, best_rp);

deflection_angle = 2*asin(1/eh);

%earth hyperbolic

a_h_earth = -mu_earth/(2.9986^2);

r_pi = a_pi*(1-e_pi);                 % Periapsis distance of the parking orbit
v_ocn_earth = sqrt(mu_earth/r_pi)*sqrt(1+e_pi);              % Orbital velocity at periapsis of the parking orbit

e_h_earth = 1 - r_pi/a_h_earth;         % Eccentricity of the hyperbolic trajectory
vph_earth = sqrt(2.9986^2 + 2*mu_earth/r_pi);

impact_n_earth = -a_h_earth*sqrt(e_h_earth^2-1);
deltaV_earth = abs(v_ocn_earth - vph_earth);


fprintf('Earth deltaV: %.4f km/s at rp = %.4f km\n', deltaV_earth, r_pi);