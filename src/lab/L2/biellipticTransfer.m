function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2] = biellipticTransfer(ai, ei, af, ef, ra_t, mu)

% Bielliptic transfer between two coplanar orbits
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

    rp_i = ai * (1 - ei); % Pericenter radius of initial orbit
    rp_f = af * (1 - ef); % Pericenter radius of final orbit

    rp_t1 = rp_i; % Pericenter radius of transfer orbit 1
    rp_t2 = rp_f; % Pericenter radius of transfer orbit 2

    a_t1 = (rp_t1 + ra_t) / 2; % Semi-major axis of transfer orbit 1
    a_t2 = (rp_t2 + ra_t) / 2; % Semi-major axis of transfer orbit 2


    vp_t1 = sqrt(mu)*sqrt((2/rp_t1) - (1/a_t1)); % Velocity at pericenter of transfer orbit 1
    vp_i = sqrt(mu)*sqrt((2/rp_i) - (1/ai)); % Velocity at pericenter of initial orbit

    va_t2 = sqrt(mu)*sqrt((2/ra_t) - (1/a_t2)); % Velocity at apocenter of transfer orbit 2
    va_t1 = sqrt(mu)*sqrt((2/ra_t) - (1/a_t1)); % Velocity at apocenter of transfer orbit 1

    vp_f = sqrt(mu)*sqrt(2/rp_f - 1/af); % Velocity at pericenter of final orbit
    vp_t2 = sqrt(mu)*sqrt(2/rp_t2 - 1/a_t2); % Velocity at pericenter of transfer orbit 2


    DeltaV1 = vp_t1 - vp_i; % First impulse for the transfer
    DeltaV2 = va_t2 - va_t1; % Second impulse for the transfer
    DeltaV3 = vp_f - vp_t2; % Third impulse for the transfer

    Deltat1 = pi*sqrt(a_t1^3/mu); % Time of flight for the first leg of the transfer
    Deltat2 = pi*sqrt(a_t2^3/mu); % Time of flight for the second leg of the transfer

end