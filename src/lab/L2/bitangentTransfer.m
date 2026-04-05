function [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)

% Bitangent transfer for elliptic orbits
%
% Syntax: [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)
%
% Inputs:
% a_i       [1x1] Initial semi-major axis                [km]
% e_i       [1x1] Initial eccentricity                   [-]
% a_f       [1x1] Final semi-major axis                  [km]
% e_f       [1x1] Final eccentricity                     [-]
% type      [char] Type of bitangent transfer
% mu        [1x1] Gravitational parameter                [km^3/s^2]
%
% Outputs:
% DeltaV1   [1x1] First impulse of the transfer          [km/s]
% DeltaV2   [1x1] Second impulse of the transfer         [km/s]
% Deltat    [1x1] Time of flight for the transfer        [s]
%
% Notes:
% - The function calculates the required impulses and time of flight for a bitangent transfer between two elliptic orbits.
% - The transfer can be of two types: type 1 (prograde) and type 2 (retrograde).
% - The function assumes impulsive maneuvers and does not account for any perturbations.
% - The orbits are assumed to be coplanar and the transfer is performed in the plane of the orbits.
%

    switch type
        case 'pa'
            rp_t = a_i*(1-e_i);                         % Perigee radius of the transfer orbit
            rp_i = rp_t;                                % Perigee radius of the initial orbit                    
            ra_t = a_f*(1+e_f);                         % Apogee radius of the transfer orbit
            ra_f = ra_t; 

            a_t = (rp_t + ra_t)/2;                      % Semi-major axis of the transfer orbit
            
            vp_t = sqrt(mu)*sqrt((2/rp_t) - (1/a_t));   % Velocity at perigee of the transfer orbit
            vp_i = sqrt(mu)*sqrt((2/rp_i) - (1/a_i));   % Velocity at perigee of the initial orbit

            va_f = sqrt(mu)*sqrt((2/ra_f) - (1/a_f));   % Velocity at apogee of the transfer orbit
            va_t = sqrt(mu)*sqrt((2/ra_t) - (1/a_t));   % Velocity at apogee of the final orbit

            DeltaV1 = abs(vp_t - vp_i);                      % First impulse for the transfer
            DeltaV2 = abs(va_f - va_t);                      % Second impulse for the transfer

        case 'ap'
            ra_t = a_i*(1+e_i);                         % Apogee radius of the transfer orbit
            ra_i = ra_t;                                % Apogee radius of the initial orbit                    
            rp_t = a_f*(1-e_f);                         % Perigee radius of the transfer orbit
            rp_f = rp_t; 

            a_t = (rp_t + ra_t)/2;                      % Semi-major axis of the transfer orbit
            
            va_t = sqrt(mu)*sqrt((2/ra_t) - (1/a_t));   % Velocity at apogee of the transfer orbit
            va_i = sqrt(mu)*sqrt((2/ra_i) - (1/a_i));   % Velocity at apogee of the initial orbit

            vp_f = sqrt(mu)*sqrt((2/rp_f) - (1/a_f));   % Velocity at perigee of the transfer orbit
            vp_t = sqrt(mu)*sqrt((2/rp_t) - (1/a_t));   % Velocity at perigee of the final orbit

            DeltaV1 = abs(va_t - va_i);                      % First impulse for the transfer
            DeltaV2 = abs(vp_f - vp_t);                      % Second impulse for the transfer

        case 'pp'
            rp_t = a_i*(1-e_i);                         % Perigee radius of the transfer orbit
            rp_i = rp_t;                                % Perigee radius of the initial orbit                    
            ra_t = a_f*(1-e_f);                         % Apogee radius of the transfer orbit
            rp_f = ra_t;                                % Second impulse for the transfer

            a_t = (rp_t + ra_t)/2;                      % Semi-major axis of the transfer orbit

            vp_t = sqrt(mu)*sqrt((2/rp_t) - (1/a_t));   % Velocity at perigee of the transfer orbit
            vp_i = sqrt(mu)*sqrt((2/rp_i) - (1/a_i));   % Velocity at perigee of the initial orbit

            vp_f = sqrt(mu)*sqrt((2/rp_f) - (1/a_f));   % Velocity at perigee of the final orbit
            va_t = sqrt(mu)*sqrt((2/ra_t) - (1/a_t));   % Velocity at apogee of the transfer orbit

            DeltaV1 = abs(vp_t - vp_i);                      % First impulse for the transfer
            DeltaV2 = abs(vp_f - va_t);                      % Second impulse for the transfer

        case 'aa'
            rp_t = a_i*(1+e_i);                         % Apogee radius of the transfer orbit
            ra_i = rp_t;                                % Apogee radius of the initial orbit
            ra_t = a_f*(1+e_f);                         % Perigee radius of the transfer orbit
            ra_f = ra_t;                                % Apogee radius of the final orbit

            a_t = (rp_t + ra_t)/2;                      % Semi-major axis of the transfer orbit

            vp_t = sqrt(mu)*sqrt((2/rp_t) - (1/a_t));   % Velocity at perigee of the transfer orbit
            va_i = sqrt(mu)*sqrt((2/ra_i) - (1/a_i));   % Velocity at apogee of the initial orbit

            va_f = sqrt(mu)*sqrt((2/ra_f) - (1/a_f));   % Velocity at apogee of the final orbit
            va_t = sqrt(mu)*sqrt((2/ra_t) - (1/a_t));   % Velocity at apogee of the transfer orbit

            DeltaV1 = abs(vp_t - va_i);                      % First impulse for the transfer
            DeltaV2 = abs(va_f - va_t);                      % Second impulse for the transfer

        otherwise
            error('Invalid transfer type. Choose from: pa, ap, pp, aa.');
    end

    Deltat = pi*sqrt(a_t^3/mu);                         % Time of flight for the transfer (half period of the transfer orbit)

end