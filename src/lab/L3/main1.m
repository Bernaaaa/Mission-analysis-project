clear; clc; close all;

% =========================================================================
% REFERENCE DATA FOR LUNAR ORBIT (a=384400, e=0.0549, mu=398600.44)
% =========================================================================
%
% CASE 1: PERIGEE TO APOGEE (th1 = 0, th2 = pi)
% Symmetry test: Must be exactly half the orbital period (T/2).
% Expected TOF: 1185921.80 seconds
% Format: 13 days, 17 hours, 25 minutes, 21.80 seconds
% 
% -------------------------------------------------------------------------
%
% CASE 2: PERIGEE TO FIRST QUADRANT (th1 = 0, th2 = pi/2)
% Velocity test: The Moon is fast near perigee. TOF is less than T/4.
% Expected TOF: 551533.26 seconds
% Format: 6 days, 9 hours, 12 minutes, 13.26 seconds
%
% -------------------------------------------------------------------------
%
% CASE 3: APOGEE TO THIRD QUADRANT (th1 = pi, th2 = 3/2*pi)
% Velocity test: The Moon is slow near apogee. TOF is more than T/4.
% Expected TOF: 634388.54 seconds
% Format: 7 days, 8 hours, 13 minutes, 8.54 seconds
%
% -------------------------------------------------------------------------
%
% CASE 4: PERIGEE CROSSING (th1 = deg2rad(350), th2 = deg2rad(10))
% Continuity test: Verifies the 2*pi jump logic across the periapsis.
% Expected TOF: 117938.01 seconds
% Format: 1 days, 8 hours, 45 minutes, 38.01 seconds
%
% =========================================================================

a = 384400; % Semi-major axis in km
e = 0.0549; % Eccentricity
th1 = deg2rad(350); % Initial true anomaly in radians
th2 = deg2rad(10); % Final true anomaly in radians
mu = 398600.44; % Standard gravitational parameter for Earth in km^3/s^2


%time1 = TOF(a, e, th1, th2, mu);
time2 = TOF_E(a, e, th1, th2, mu);
time3 = TOF_M(a, e, th1, th2, mu);

total_time = sqrt(a^3 / mu) * 2 * pi; % Total orbital period
fprintf('\nTotal Orbital Period: %.2f seconds\n', total_time);
fprintf('Total Orbital Period: %d days, %d hours, %d minutes, %.2f seconds\n', floor(total_time / 86400), floor(mod(total_time, 86400) / 3600), floor(mod(total_time, 3600) / 60), mod(total_time, 60));
