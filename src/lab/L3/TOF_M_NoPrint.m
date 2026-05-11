function deltat = TOF_M_NoPrint(a, e, th1, th2, mu)
% Time of Flight 
%
% Syntax: deltat = TOF_M_NoPrint(a, e, th1, th2, mu)
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
% - The time of flight is computed using the difference in mean anomalies corresponding to the two true anomalies, divided by the mean motion of the orbit.
    

    if e < 0 || e >= 1
        error('The only accepted values of eccentricity are in the range [0, 1) for elliptical orbits.');
    end


    E1 = 2 * atan2(sqrt(1-e) * sin(th1/2), sqrt(1+e)*cos(th1/2));       % Eccentric anomaly corresponding to th1
    E2 = 2 * atan2(sqrt(1-e) * sin(th2/2), sqrt(1+e)*cos(th2/2));       % Eccentric anomaly corresponding to th2

    M1 = E1 - e*sin(E1);    % Mean anomaly corresponding to th1
    M2 = E2 - e*sin(E2);    % Mean anomaly corresponding to th2


    if M2 < M1              % Ensure that the time of flight is positive (add 2*pi if th2 crosses perigee in the next orbit)
        M2 = M2 + 2*pi;
    end

    dM = M2 - M1;           % Difference in mean anomalies
    n = sqrt(mu / a^3);     % Mean motion
    deltat =  dM / n;       % Time of flight
end