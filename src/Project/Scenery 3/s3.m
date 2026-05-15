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
M_sun = 1.989e+30;                      % Mass of the Sun [kg]
AU = 149597870.7;                       % Astronomical Unit [km]
M = 1.0472E+12;                         % Mass of 3908 Nyx [kg]
D = 1.00;                               % Diameter of 3908 Nyx [km]
G = 6.67430e-20;                        % Gravitational constant [km^3/kg/s^2]
mu = G*M;                               % Gravitational parameter of 3908 Nyx [km^3/s^2]
    % Sphere of Influence of 3908 Nyx [km] %CHECK FOR VALIDITY

%% 2. ORBITAL ELEMENTS OF 3908 NYX
a = 1.927925;
e_nyx = 0.459072;
i_nyx = deg2rad(2.19);
OM_nyx = deg2rad(261.19);
om_nyx = deg2rad(126.66);
M_nyx = 1.0472E+12;
D_nyx = 1.00;

a_nyx  = a * AU;                 % Semi-major axis [km]
v_escape_nyx = sqrt(2*mu/(D_nyx/2));     % Escape velocity from 3908 Nyx [km/s]



%% 3 Earth's orbital elements

a_ea  = 1.4946e+08;                  % Semi-major axis [km]
e_ea  = 0.016;                       % Eccentricity [-]
i_ea  = 9.1920e-05;                  % Inclination [rad]
OM_ea = 2.7847;                      % Longitude of Ascending Node [rad]
om_ea = 5.2643;                      % Argument of Pericenter [rad]

mu_ea = 398600.4418;             % Gravitational parameter of Earth [km^3/s^2]
r_ea = 6371;                      % Radius of Earth [km]


%% ORBITAL ELEMENTS OF PARKING EARTH ORBIT
r_xf = -6103.007500;
r_yf = 13604.661000;
r_zf = 7991.352600;

v_xf = -5.163000;
v_yf = -3.151000;
v_zf = 1.421000;


rf = [r_xf; r_yf; r_zf];
vf = [v_xf; v_yf; v_zf];


[a_start, e_start, i_start, OM_start, om_start, th_start] = car2par(rf, vf, mu);


%% FIRST HYPERBOLIC ESCAPE FROM EARTH
% Vinf = %Vtransfer s2 - V earth (prendere da output s2)


