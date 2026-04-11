function deltat = TOF_E(a, e, th1, th2, mu)
% Time of Flight 
%
% Syntax: deltat = TOF_E(a, e, th1, th2, mu)
%
% Inputs:
% a           [1x1]    Semi-major axis                             [km]
% e           [1x1]    Eccentricity                                [-]
% th1         [1x1]    Initial true anomaly                        [rad]
% th2         [1x1]    Final true anomaly                          [rad]
% mu          [1x1]    Gravitational parameter                     [km^3/s^2]
%
% Outputs:
% deltat      [1x1]    Time of flight between th1 and th2          [s]
%
% Notes: 
% - The function calculates the time of flight (dt) for a spacecraft to travel between two true anomalies (phi_1 and phi_2) in an elliptical orbit defined by its semi-major axis (a) and eccentricity (e).
% - The time of flight is computed using the difference in eccentric anomalies, modified to ensure a forward temporal direction.

    % Check for orbit validity (elliptical only)
    if e < 0 || e >= 1
        error('The only accepted values of eccentricity are in the range [0, 1) for elliptical orbits.');
    end

    % Compute Eccentric Anomalies using the half-angle tangent identity.
    % atan2 is used to handle the correct quadrant based on the sign of sin/cos.
    E1 = 2 * atan2(sqrt(1-e) * sin(th1/2), sqrt(1+e) * cos(th1/2));
    E2 = 2 * atan2(sqrt(1-e) * sin(th2/2), sqrt(1+e) * cos(th2/2));

    % Ensure that the time of flight is positive.
    if E2 < E1
        E2 = E2 + 2 * pi;
    end

    % Calculate the time of flight using Kepler's Equation expanded form.
    deltat = sqrt(a^3 / mu) * ((E2 - E1) - e * (sin(E2) - sin(E1)));

    % Display the time of flight in seconds and also in (days, hours, minutes, seconds) format.
    fprintf('\nTime of Flight: %.2f seconds\n', deltat);

    days = floor(deltat / 86400);
    hours = floor(mod(deltat, 86400) / 3600);
    minutes = floor(mod(deltat, 3600) / 60);
    seconds = mod(deltat, 60);

    fprintf('Time of Flight: %d days, %d hours, %d minutes, %.2f seconds\n\n', days, hours, minutes, seconds);

end