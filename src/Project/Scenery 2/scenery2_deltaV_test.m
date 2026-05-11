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

[best_th_earth, best_th_nyx, best_om_tn, best_dv] = grid_search_s2(a_i, e_i, i_i, OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu, R_sun, grid_size);

x0 = [best_th_earth, best_th_nyx, best_om_tn];
lb = [0, 0, 0];
ub = [2*pi, 2*pi, 2*pi];

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ... 
    'TolFun', 1e-10, ...
    'MaxIterations', 500);

[x_opt, dv_opt] = fmincon(@(x) mission_transfer_nyx(x, a_i, e_i, i_i, OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu, R_sun), ...
    x0, [], [], [], [], lb, ub, [], options);

% Stampiamo il confronto
fprintf('\n==========================================\n');
fprintf('RISULTATO GRID SEARCH: %.6f km/s\n', best_dv);
fprintf('RISULTATO FMINCON:     %.6f km/s\n', dv_opt);
fprintf('Guadagno ottenuto:     %.6f km/s\n', best_dv - dv_opt);
fprintf('==========================================\n');

th_earth_opt = x_opt(1);
th_nyx_opt   = x_opt(2);
om_tn_opt    = x_opt(3);

% 2. Ricalcoliamo i vettori posizione e velocità ottimi
[r1_opt, v_earth_opt] = par2car(a_i, e_i, i_i, OM_i, om_i, th_earth_opt, mu);
[r2_opt, v_nyx_opt]   = par2car(a_f, e_f, i_f, OM_f, om_f, th_nyx_opt, mu);

% 3. Rieseguiamo la geometria del piano per questi vettori
p_vect_opt = cross(r1_opt, r2_opt);
h_t_opt = p_vect_opt / norm(p_vect_opt);

% Inclinazione e Nodo Ascendente dell'orbita ottima
i_t_opt = acos(h_t_opt(3));
p_kht_opt = cross([0; 0; 1], p_vect_opt);
N_t_opt = p_kht_opt / norm(p_kht_opt);
if N_t_opt(2) >= 0
    OM_t_opt = acos(N_t_opt(1));
else
    OM_t_opt = 2*pi - acos(N_t_opt(1));
end

% Calcolo delle anomalie vere sul trasferimento (th_i_t e th_f_t)
% Dobbiamo rifare la rotazione sul piano per trovare alpha1 e alpha2
R3_OM_opt = [cos(OM_t_opt) sin(OM_t_opt) 0; -sin(OM_t_opt) cos(OM_t_opt) 0; 0 0 1];
R1_i_opt  = [1 0 0; 0 cos(i_t_opt) sin(i_t_opt); 0 -sin(i_t_opt) cos(i_t_opt)];
r1_plane_opt = (R1_i_opt * R3_OM_opt) * r1_opt;
r2_plane_opt = (R1_i_opt * R3_OM_opt) * r2_opt;

alpha1_opt = atan2(r1_plane_opt(2), r1_plane_opt(1));
alpha2_opt = atan2(r2_plane_opt(2), r2_plane_opt(1));
while alpha2_opt < alpha1_opt, alpha2_opt = alpha2_opt + 2*pi; end

th_i_t_opt = alpha1_opt - om_tn_opt;
th_f_t_opt = alpha2_opt - om_tn_opt;

% 4. Parametri orbitali finali
r1_n = norm(r1_opt);
r2_n = norm(r2_opt);
e_t_opt = (r2_n - r1_n) / (r1_n*cos(th_i_t_opt) - r2_n*cos(th_f_t_opt));
p_t_opt = r1_n * (1 + e_t_opt*cos(th_i_t_opt));
a_t_opt = p_t_opt / (1 - e_t_opt^2);

% 5. STAMPA DEI RISULTATI FINALI
fprintf('\n--- PARAMETRI ORBITALI OTTIMIZZATI ---\n');
fprintf('Semiasse Maggiore:  %.2f km\n', a_t_opt);
fprintf('Eccentricità:       %.6f\n', e_t_opt);
fprintf('Inclinazione:       %.4f gradi\n', rad2deg(i_t_opt));
fprintf('Long. Nodo Asc.:    %.4f gradi\n', rad2deg(OM_t_opt));
fprintf('Arg. Pericentro:    %.4f gradi\n', rad2deg(om_tn_opt));
fprintf('Anomalia Vera Partenza: %.4f gradi\n', rad2deg(th_i_t_opt));
fprintf('Anomalia Vera Arrivo:   %.4f gradi\n', rad2deg(th_f_t_opt));


% --- CALCOLO VETTORI VELOCITÀ OTTIMIZZATI PER LA STAMPA ---

% Calcolo velocità nel piano perifocale (usando i dati opt calcolati prima)
v_coeff_opt = sqrt(mu/p_t_opt);
v1t_pf_opt = v_coeff_opt * [-sin(th_i_t_opt); e_t_opt + cos(th_i_t_opt); 0];
v2t_pf_opt = v_coeff_opt * [-sin(th_f_t_opt); e_t_opt + cos(th_f_t_opt); 0];

% Matrici di rotazione finale (per tornare in Eliocentrico)
R3_om_opt_T = [cos(om_tn_opt) -sin(om_tn_opt) 0; sin(om_tn_opt) cos(om_tn_opt) 0; 0 0 1];
RT_Plane_to_ELIO_opt = (R1_i_opt * R3_OM_opt)'; % Inversa della rotazione al piano

% Vettori velocità finali nel sistema Eliocentrico
v1t_final = RT_Plane_to_ELIO_opt * (R3_om_opt_T * v1t_pf_opt);
v2t_final = RT_Plane_to_ELIO_opt * (R3_om_opt_T * v2t_pf_opt);

% --- STAMPA FINALE DEI VETTORI ---
fprintf('\n--- VETTORI VELOCITÀ OTTIMIZZATI (km/s) ---\n');
fprintf('Velocità Terra alla partenza (v_i):   [%.4f, %.4f, %.4f]\n', v_earth_opt);
fprintf('Velocità Trasferimento partenza (v1t): [%.4f, %.4f, %.4f]\n', v1t_final);
fprintf('DeltaV1 (Partenza):                   %.4f km/s\n\n', norm(v1t_final - v_earth_opt));

fprintf('Velocità Nyx all''arrivo (v_f):        [%.4f, %.4f, %.4f]\n', v_nyx_opt);
fprintf('Velocità Trasferimento arrivo (v2t):   [%.4f, %.4f, %.4f]\n', v2t_final);
fprintf('DeltaV2 (Arrivo):                     %.4f km/s\n', norm(v_nyx_opt - v2t_final));
fprintf('------------------------------------------\n');
T = TOF_M(a_t_opt, e_t_opt, th_i_t_opt, th_f_t_opt, mu);
fprintf('Time of Flight: %.2f seconds (%.2f days)\n', T, T/86400);