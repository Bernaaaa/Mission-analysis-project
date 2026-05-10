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

grid_size = 100;

th_i = linspace(0, 2*pi, grid_size); % true anomaly of the initial orbit varying from 0 to 360 degrees
th_f = linspace(0, 2*pi, grid_size); % true anomaly of the final orbit varying from 0 to 360 degrees
om_tn = linspace(0, 2*pi, grid_size); % pericenter argument of the transfer orbit varying from 0 to 360 degrees

r1 = zeros(3, grid_size);
v_i = zeros(3, grid_size);
r2 = zeros(3, grid_size);
v_f = zeros(3, grid_size);

r1_n = zeros(1, grid_size);
r2_n = zeros(1, grid_size);

best_dv = inf;

for i = 1:grid_size
    [r1(:,i), v_i(:,i)] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i(i), mu);
    [r2(:,i), v_f(:,i)] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f(i), mu); 
    r1_n(i) = norm(r1(:,i));
    r2_n(i) = norm(r2(:,i));
end

for j = 1:grid_size

    r1_current = r1(:,j);  v_i_current = v_i(:,j);  r1_n_current = r1_n(j);

    for k = 1:grid_size

        r2_current = r2(:,k);  v_f_current = v_f(:,k);  r2_n_current = r2_n(k);
                        
        p_vect_r1r2 = cross(r1_current, r2_current);
        h_t = p_vect_r1r2 / norm(p_vect_r1r2);
        i_t = acos(h_t(3));

        p_vect_kht = cross([0; 0; 1], p_vect_r1r2);
        N_t = p_vect_kht / norm(p_vect_kht);

        if N_t(2) >= 0
            OM_t = acos(N_t(1));
        else
            OM_t = 2*pi - acos(N_t(1));
        end

        R3_OM = [ cos(OM_t)  sin(OM_t)  0;
                 -sin(OM_t)  cos(OM_t)  0;
                     0          0       1];

        R1_i  = [ 1          0          0;
                  0       cos(i_t)   sin(i_t);
                  0      -sin(i_t)   cos(i_t)];

        RT_ELIO_to_Plane = R1_i * R3_OM;

        r1_plane = RT_ELIO_to_Plane * r1_current;
        r2_plane = RT_ELIO_to_Plane * r2_current;

        alpha1 = atan2(r1_plane(2), r1_plane(1));
        alpha2 = atan2(r2_plane(2), r2_plane(1));

        RT_Plane_to_ELIO = RT_ELIO_to_Plane';

        for om_tn_val = om_tn

            th_i_tn = alpha1 - om_tn_val;
            th_f_tn = alpha2 - om_tn_val;

            num_e = r2_n_current - r1_n_current;
            den_e = r1_n_current*cos(th_i_tn) - r2_n_current*cos(th_f_tn);
            e_t = num_e / den_e;

            if e_t < 0 || e_t >= 1, continue; end % Salto ellissi non valide

            num_a = r1_n_current * (1 + e_t*cos(th_i_tn));
            den_a = 1 - e_t^2;
            a_t = num_a / den_a;           
            rp_t = a_t * (1 - e_t);

            if rp_t <= R_sun + 1000, continue; end

            % Calcolo velocità (Logica par2car integrata ed efficiente)
            p_t = a_t * (1 - e_t^2);
            v_coeff = sqrt(mu/p_t);

            % Velocità nel frame perifocale (PF)
            v1t_pf = v_coeff * [-sin(th_i_tn); e_t + cos(th_i_tn); 0];
            v2t_pf = v_coeff * [-sin(th_f_tn); e_t + cos(th_f_tn); 0];

            R3_om = [ cos(om_tn_val)  sin(om_tn_val)  0;
                     -sin(om_tn_val)  cos(om_tn_val)  0;
                            0               0         1];

            R3_om_T = R3_om';

            v1t = RT_Plane_to_ELIO * (R3_om_T * v1t_pf);
            v2t = RT_Plane_to_ELIO * (R3_om_T * v2t_pf);     

            DeltaV = norm(v1t - v_i_current) + norm(v_f_current - v2t);
                
            if DeltaV < best_dv
                best_dv = DeltaV;
                best_v1t = v1t;
                best_v2t = v2t;
                best_v_i = v_i_current;
                best_v_f = v_f_current;
                best_a_t = a_t;
                best_e_t = e_t;
                best_i_t = i_t;
                best_OM_t = OM_t;
                best_om_tn = om_tn_val;
                best_th_i_tn = th_i_tn;
                best_th_f_tn = th_f_tn;
            end
        end
    end
end

% Print starting and ending velocities, and best transfer orbit parameters
fprintf('Initial velocity (v_i): [%.4f, %.4f, %.4f] km/s\n', best_v_i);
fprintf('Final velocity (v_f): [%.4f, %.4f, %.4f] km/s\n', best_v_f);
fprintf('Transfer velocity at departure (v1t): [%.4f, %.4f, %.4f] km/s\n', best_v1t);
fprintf('Transfer velocity at arrival (v2t): [%.4f, %.4f, %.4f] km/s\n\n', best_v2t);


fprintf('Best DeltaV: %.4f km/s\n\n', best_dv);

% Detailed transfer orbit parameters
fprintf('Best Transfer Orbit Parameters:\n');
fprintf('Semi-major axis: %.4f km\n', best_a_t);
fprintf('Eccentricity: %.4f\n', best_e_t);
fprintf('Inclination: %.4f degrees\n', rad2deg(best_i_t));
fprintf('Longitude of ascending node: %.4f degrees\n', rad2deg(best_OM_t));
fprintf('Argument of periapsis: %.4f degrees\n', rad2deg(best_om_tn));
fprintf('True anomaly at initial position: %.4f degrees\n', rad2deg(best_th_i_tn));
fprintf('True anomaly at final position: %.4f degrees\n', rad2deg(best_th_f_tn));
