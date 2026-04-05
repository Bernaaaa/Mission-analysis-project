function [DeltaV, thi, thf] = changePericenterArg(a, e, omi, omf, mu)

% changePericenterArg - Change the argument of pericenter of an orbit
%
% Syntax: [DeltaV, thi, thf] = changePericenterArg(a, e, omi, omf, mu)
%
% Inputs:
% a         [1x1] Semi-major axis of the orbit             [km]
% e         [1x1] Eccentricity of the orbit                [-]
% omi       [1x1] Initial argument of pericenter           [rad]
% omf       [1x1] Final argument of pericenter             [rad]
% mu        [1x1] Gravitational parameter                  [km^3/s^2]
%
% Outputs:
% DeltaV    [1x1] Maneuver impulse                         [km/s]
% thi       [2x1] Initial true anomalies of the maneuver     [rad]
% thf       [2x1] Final true anomalies of the maneuver       [rad]
%
% Notes:
% - The maneuver is performed at the point where the argument of pericenter is changed.
% - The function assumes impulsive maneuvers and does not account for any perturbations.
% - The orbit's shape and size remain unchanged, only the orientation of the pericenter is altered.
% - The change can be performed at two true anomalies.
%

    thi = [0, 0];
    thf = [0, 0];

    deltaO = omf - omi;                             % Change in argument of pericenter

    thi(1) = deltaO/2; 
    thi(2) = thi(1) + pi;

    thf(1) = 2*pi - thi(1);
    thf(2) = pi - thi(1);

    p = a*(1-e^2);                                  % Semi-latus rectum

    DeltaV = 2*sqrt(mu/p)*e*sin(deltaO/2);          % Impulse required for the maneuver

end