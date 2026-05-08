clear; clc; close all;

%: Heliocentric transfer orbit characterization from earth to asteroid 3908 Nyx
% a_f is in astronomical units, but we need it in kilometers for the calculations, so we will convert it using the conversion factor 1 AU = 149597870.7 km.

% Asteroid 3908 Nyx (1980 PA) data from JPL Small-Body Database:
a_f = 1.927925 * 149597870.7; % km
e_f = 0.459072; 
i_f = deg2rad(2.19);
OM_f = deg2rad(261.19);     
om_f = deg2rad(126.66);  
M = 1.0472E+12; 
D = 1.00;

% Initial orbit parameters (Earth's orbit around the Sun):
a_i = 1.4946e+08; % km
e_i = 0.016;
i_i = 9.1920e-05; % radians
OM_i = 2.7847; % radians
om_i = 5.2643; % radians

mu = 1.32712440041279419e+11; % km^3/s^2 from https://ssd.jpl.nasa.gov/astro_par.html (NASA JPL)
R_sun = 696340; % km

th_i = linspace(0, 2*pi, 100); % true anomaly of the initial orbit varying from 0 to 360 degrees
th_f = linspace(0, 2*pi, 100); % true anomaly of the final orbit varying from 0 to 360 degrees
om_tn = linspace(0, 2*pi, 100); % pericenter argument of the transfer orbit varying from 0 to 360 degrees
best_dv = inf;

for om_tn_val = om_tn
    for th_i_val = th_i
        for th_f_val = th_f
            try 
                [r1_i, v_i] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i_val, mu);
                [r2_f, v_f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f_val, mu);

                r1 = r1_i;
                r2 = r2_f;
                
                p_vect_r1r2 = cross(r1, r2);
                p_vect_kht = cross([0; 0; 1], p_vect_r1r2);

                h_t = p_vect_r1r2 / norm(p_vect_r1r2);
                i_t = acos(h_t(3));
                N_t = p_vect_kht / norm(p_vect_kht);

                if N_t(2) >= 0
                    OM_t = acos(N_t(1));
                else
                    OM_t = 2*pi - acos(N_t(1));
                end

                % Transfer matrix from 𝑇𝐸𝑙𝑖𝑜→𝑃𝐹𝑇𝑛 = 𝑅𝜔𝑇𝑛 𝑅𝑖𝑇 𝑅Ω𝑇
                R3_OM = [ cos(OM_t)  sin(OM_t)  0;
                        -sin(OM_t)  cos(OM_t)  0;
                        0          0          1];

                R1_i  = [ 1          0          0;
                        0          cos(i_t)   sin(i_t);
                        0         -sin(i_t)   cos(i_t)];

                R3_om = [ cos(om_tn_val)  sin(om_tn_val)  0;
                        -sin(om_tn_val)  cos(om_tn_val)  0;
                        0               0               1];

                RT_EL_PF = R3_om * R1_i * R3_OM;

                r1_tn = RT_EL_PF * r1;
                r2_tn = RT_EL_PF * r2;

                
                cos_th_i = r1_tn(1) / norm(r1_tn);
                sin_th_i = r1_tn(2) / norm(r1_tn);

                cos_th_f = r2_tn(1) / norm(r2_tn);
                sin_th_f = r2_tn(2) / norm(r2_tn);

                % Calculate the true anomalies of the initial and final orbits in the transfer plane
                th_i_tn = atan2(sin_th_i, cos_th_i);
                th_f_tn = atan2(sin_th_f, cos_th_f);

                
                num_e = norm(r2) - norm(r1);
                den_e = norm(r1)*cos(th_i_tn) - norm(r2)*cos(th_f_tn);
                e_t = num_e / den_e;

                num_a = norm(r1) * (1 + e_t*cos(th_i_tn));
                den_a = 1 - e_t^2;
                a_t = num_a / den_a;

                rp_t = a_t * (1 - e_t);

                if ~(0 <= e_t && e_t < 1) || ~(rp_t > R_sun + 1000)
                    continue; % Skip if the eccentricity or periapsis is not valid
                end

                h_earth = cross(r1, v_i);
                

                [~, v1t] = par2car(a_t, e_t, i_t, OM_t, om_tn_val, th_i_tn, mu);
                [~, v2t] = par2car(a_t, e_t, i_t, OM_t, om_tn_val, th_f_tn, mu);

                DeltaV = abs(norm(v1t - v_i)) + abs(norm(v_f - v2t));
                
                if DeltaV < best_dv
                    best_dv = DeltaV;
                    best_a_t = a_t;
                    best_e_t = e_t;
                    best_i_t = i_t;
                    best_OM_t = OM_t;
                    best_om_tn = om_tn_val;
                    best_th_i_tn = th_i_tn;
                    best_th_f_tn = th_f_tn;
                end
            catch
                continue; % If the geometry is impossible, skip to the next iteration
            end
        end
    end
end

fprintf('Best DeltaV: %.4f km/s\n', best_dv);
fprintf('Best Transfer Orbit Parameters:\n');
fprintf('Semi-major axis: %.4f km\n', best_a_t);
fprintf('Eccentricity: %.4f\n', best_e_t);
fprintf('Inclination: %.4f degrees\n', rad2deg(best_i_t));
fprintf('Longitude of ascending node: %.4f degrees\n', rad2deg(best_OM_t));
fprintf('Argument of periapsis: %.4f degrees\n', rad2deg(best_om_tn));
fprintf('True anomaly at initial position: %.4f degrees\n', rad2deg(best_th_i_tn));
fprintf('True anomaly at final position: %.4f degrees\n', rad2deg(best_th_f_tn));
