function [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)
% Change plane of maneuver
%
% Syntax: [DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)
%
% Inputs:
% a           [1x1]    Semi-major axis                              [km]
% e           [1x1]    Eccentricity                                 [-]
% i_i         [1x1]    Initial inclination                          [rad]
% OMi         [1x1]    Initial right ascension of ascending node    [rad]
% omi         [1x1]    Initial pericenter anomaly                   [rad]
% i_f         [1x1]    Final inclination                            [rad]
% OMf         [1x1]    Final right ascension of ascending node      [rad]
% mu          [1x1]    Gravitational parameter                      [km^3/s^2]
%
% Outputs:
% DeltaV       [1x1]    Required velocity change for the maneuver   [km/s]
% omf          [1x1]    Final pericenter anomaly                    [rad]
% theta        [1x1]    True anomaly at maneuver                    [rad]
%
% Notes: - The function calculates the required DeltaV to change the orbital plane from an initial configuration (i_i, OMi) 
% to a final configuration (i_f, OMf) while keeping the same semi-major axis and eccentricity. 
% It also determines the true anomaly (theta) at which the maneuver should be performed and the 
% new argument of pericenter (omf) after the maneuver.
% 

    dOM = OMf - OMi;
    dI  = i_f - i_i;
    
    cos_alpha = cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(dOM);    % Law of cosines for spherical triangles
    alpha = acos(max(-1, min(1, cos_alpha)));   % Numerical safety for acos
    if alpha < 1e-12
        DeltaV = 0; omf = omi; theta = 0; return;
    end

    sin_ui = (sin(dOM) * sin(i_f)) / sin(alpha); % Law of sines for spherical triangles
    sin_uf = (sin(dOM) * sin(i_i)) / sin(alpha); % Law of sines for spherical triangles

    cos_ui = (cos(i_f) - cos(alpha)*cos(i_i)) / (sin(alpha)*sin(i_i)); % Law of cosines for spherical triangles
    cos_uf = (cos(i_i) - cos(alpha)*cos(i_f)) / (sin(alpha)*sin(i_f)); % Law of cosines for spherical triangles

    if (dOM * dI) > 0                           % Same direction of change in RAAN and inclination
        ui = atan2(sin_ui, -cos_ui);            % Argument of latitude of maneuver point in initial orbit
        uf = atan2(sin_uf,  cos_uf);            % Argument of latitude of maneuver point in final orbit
        theta = ui - omi;                       % True anomaly of maneuver point in initial orbit
        omf = uf - theta;                       % Argument of pericenter in final orbit
    else
        ui = atan2(sin_ui,  cos_ui);
        uf = atan2(sin_uf, -cos_uf);
        theta = 2*pi - ui - omi;
        omf = 2*pi - uf - theta;
    end

    theta = mod(theta, 2*pi);                   % Ensure theta is in [0, 2*pi]
    omf   = mod(omf, 2*pi);                     % Ensure omf is in [0, 2*pi]


    p = a * (1 - e^2);
    r1 = p / (1 + e * cos(theta));    
    v_theta = sqrt(mu/p) * (1 + e * cos(theta));
    DeltaV = 2 * v_theta * sin(alpha / 2);

    fprintf('Calculated Maneuver:\n');
    fprintf('Theta (Manuever Point): %.2f deg\n', rad2deg(theta));
    fprintf('New Argument of Pericenter (omf): %.2f deg\n', rad2deg(omf));
    fprintf('DeltaV Required: %.4f km/s\n\n', DeltaV);

    r2 = p / (1 + e * cos(theta + pi));  

    if r2 > r1 + 1e-6                           % If the opposite point is farther, check if it offers a better maneuver

        theta = mod(theta + pi, 2*pi);
        v_theta = sqrt(mu/p) * (1 + e * cos(theta));
        DeltaV1 = 2 * v_theta * sin(alpha / 2);

        fprintf('Since the opposing point is farther away, it is more convenient to maneuver at:\ntheta* = theta + pi = %.2f deg\n', rad2deg(theta));
        fprintf('Maneuver cost at theta* = %.4f km/s\n', DeltaV1);
        fprintf('Earned DeltaV = %.4f km/s\n\n', DeltaV - DeltaV1);

        DeltaV = DeltaV1;
    end

end