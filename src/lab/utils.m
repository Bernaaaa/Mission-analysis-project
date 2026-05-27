% Costanti e parametri iniziali
mu_earth = 398600.4418;
AU = 149597870.7;                   % Astronomical Unit [km]
mu_sun = 1.32712440041279419e+11;    % Solar Gravitational Parameter [km^3/s^2]
R_sun = 696340;                     % Sun Radius [km]
day2sec = 86400;                    % Conversion factor: days to seconds
M_sun   = 1.989e30;       % [kg] Mass of the Sun
M_earth = 5.972e24;       % [kg] Mass of the Earth
M_nyx   = 1.0472e12;      % [kg] Mass of Asteroid 3908 Nyx
D_nyx   = 1.0;            % [km] Diameter of 3908 Nyx
R_earth = 6378.137;       % [km] Earth's equatorial radius
G        = 6.67430e-20;   
mu_nyx   = G * M_nyx;     % [km^3/s^2] Nyx's gravitational parameter

% Parametri orbitali iniziali (GTO)
GTO.a = 24400.00;
GTO.e = 0.728300;
GTO.i = 0.104700;
GTO.OM = 2.361000;
GTO.om = 3.107000;
GTO.th = 2.135000;
GTO.mu = mu_earth;

% Posizione e velocità finali (orbita di parcheggio)
Park.xf = -7090.590200;
Park.yf = -5612.557300;
Park.zf = 3948.902900;
Park.vx = 5.698000;
Park.vy = -5.995000;
Park.vz = 1.710000;

% Parametri orbitali della terra
Earth.a = 1.4946e+08;                  % Semi-major axis [km]
Earth.e = 0.016;                       % Eccentricity [-]
Earth.i = 9.1920e-05;                  % Inclination [rad]
Earth.OM = 2.7847;                      % Longitude of Ascending Node [rad]
Earth.om = 5.2643;  

% Parametri orbitali del target (asteroide 3908 Nyx)
Nyx.a = 1.927925 * AU;               % Semi-major axis [km]
Nyx.e = 0.459072;                    % Eccentricity [-]
Nyx.i = deg2rad(2.19);               % Inclination [rad]
Nyx.OM = deg2rad(261.19);             % Longitude of Ascending Node [rad]
Nyx.om = deg2rad(126.66);             % Argument of Perihelion [rad]

