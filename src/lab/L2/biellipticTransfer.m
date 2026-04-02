function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(ai, ei, af, ef, ra_t, mu)

% Bitangent transfer between two coplanar orbits
%
% Inputs:
% ai            [1x1]   Initial orbit semi-major axis       [km]
% ei            [1x1]   Initial orbit eccentricity          [-]
% af            [1x1]   Final orbit semi-major axis         [km]
% ef            [1x1]   Final orbit eccentricity            [-]
% ra_t          [1x1]   Transfer orbit apocenter distance   [km]
% mu            [1x1]   Gravitational parameter             [km^3/s^2]
%
% Outputs:
% DeltaV1       [1x1]   First impulse for the transfer                      [km/s]
% DeltaV2       [1x1]   Second impulse for the transfer                     [km/s]
% DeltaV3       [1x1]   Third impulse for the transfer                      [km/s]
% Deltat1       [1x1]   Time of flight for the first leg of the transfer    [s]
% Deltat2       [1x1]   Time of flight for the second leg of the transfer   [s]
%




    
end