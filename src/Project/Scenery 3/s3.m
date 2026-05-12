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
a_nyx  = 1.927925 * AU;                 % Semi-major axis [km]

v_escape = sqrt(2*mu/(D/2));            % Escape velocity from 3908 Nyx [km/s]

r_SOI_nyx = a_nyx * (M/M_sun)^(2/5);    % Sphere of Influence of 3908 Nyx [km] %CHECK FOR VALIDITY

