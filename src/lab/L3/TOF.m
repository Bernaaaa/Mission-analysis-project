function deltat = TOF(a, e, th1, th2, mu)
% Time of Flight 
%
%deltat = TOF(a, e, th1, th2, mu)
%
% -----------------------------------------------------------------
%Inputs arguments:
% a           [1x1]    Semi-major axis                             [km]
% e           [1x1]    Eccentricity                                [-]
% th1         [1x1]    Initial true anomaly                        [rad]
% th2         [1x1]    Final true anomaly                          [rad]
% mu          [1x1]    Gravitational parameter                     [km^3/s^2]
%
%-------------------------------------------------------------------
%Outputs arguments:
% deltat      [1x1]    Time of flight between th1 and th2
%
%
% Notes: - The function calculates the time of flight (dt) for a spacecraft to travel between two true anomalies (phi_1 and phi_2) in an elliptical orbit defined by its semi-major axis (a) and eccentricity (e).
%        - The time of flight is computed using the difference in mean anomalies corresponding to the two true anomalies, divided by the mean motion of the orbit.


if e >=1 || e <= 0
    error('Eccentricity must be in the range [0, 1) for elliptical orbits.');
end

E1 = 2 * atan2( sqrt((1-e)/(1+e)) * tan(th1 /2) );   % Eccentric anomaly at th1
E2 = 2 * atan2( sqrt((1-e)/(1+e)) * tan(th2 /2) );   % Eccentric anomaly at th2

deltat = sqrt(a^3 / mu) * (E2 - E1 - e * (sin(E2) - sin(E1)));  % Time of flight

if th2 < th1   
    T = 2 * pi * sqrt(a^3 / mu);  % Orbital period  
    deltat = deltat + T;  % Adjust for full orbit if th2 is less than th1
end



