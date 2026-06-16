% =========================================================================
% INTERPLANETARY MISSION ANALYSIS - SCENARIO 3 (AVANZATO)
% Calcolo dell'Iperbole di Fuga non complanare
% =========================================================================

clear; clc; close all;

%% 0. PATH SETUP E COSTANTI FISICHE
currentDir = fileparts(mfilename('fullpath'));
funcFolder = fullfile(currentDir, '..', '..', 'lab');
if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
end

% Costanti
mu = 398600.4418;       % [km^3/s^2] Parametro gravitazionale Terra
R_earth = 6378.1363;    % [km] Raggio medio terrestre
epsilon = deg2rad(23.45); % [rad] Obliquità dell'eclittica

opts_fsolve = optimoptions('fsolve', 'Display', 'off', ...
    'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);

%% 1. ORBITA DI PARCHEGGIO (Da Scenario 1)
% Vettore di stato ECI iniziale
ri = [-7090.590200; -5612.557300;  3948.902900];
vi = [   5.698000;    -5.995000;     1.710000];

[a_i, e_i, i_i, OM_i, om_i, theta_iniziale] = car2par(ri, vi, mu);

%% 2. VELOCITÀ DI ECCESSO IPERBOLICO (Da Scenario 2)
% Dati in sistema Eclittico
V_Earth_Depart =  [ 2.9766,  29.5101, -0.0026];
V_Trans_Depart =  [ 3.1273,  32.4459,  0.5897];

% V_inf in Eclittica
v_inf_earth = V_Trans_Depart - V_Earth_Depart;


T_ECI_to_HELIO = [1, 0,            0; 
             0, cos(epsilon), sin(epsilon); 
             0, -sin(epsilon), cos(epsilon)];

v_inf_earth_eci = (T_ECI_to_HELIO' * v_inf_earth')';
v_inf_norm = norm(v_inf_earth_eci);

% Versore asintotico di uscita (r_inf)
r_inf = (v_inf_earth_eci / v_inf_norm)';

% Semiasse maggiore dell'iperbole (dalla conservazione dell'energia)
a_h = -mu / (v_inf_norm^2);


%% 3. OTTIMIZZAZIONE DEL PUNTO DI INIEZIONE (Grid Search)
fprintf('Inizio ottimizzazione del punto di iniezione (Grid Search su 360 gradi)...\n');

theta_val = 0 : deg2rad(1) : 2*pi;
DV_values = inf(size(theta_val));
exit_flags = zeros(size(theta_val));

e_h_vals  = zeros(size(theta_val));
i_H_vals  = zeros(size(theta_val));
OM_H_vals = zeros(size(theta_val));
om_H_vals = zeros(size(theta_val));

for j = 1:length(theta_val)
    
    th_in = theta_val(j);
    
    % Posizione e velocità sull'orbita di parcheggio
    [rh_vec, v_park] = par2car(a_i, e_i, i_i, OM_i, om_i, th_in, mu);
    r_h_norm = norm(rh_vec);
    rh_ver = rh_vec / r_h_norm;
    
    % Angolo di separazione
    alpha = acos(dot(rh_ver, r_inf));
    
    % Risoluzione per l'eccentricità
    g = @(e) (a_h * (1 - e.^2)) ./ (1 + e .* cos(acos(-1 ./ e) - alpha)) - r_h_norm;
    [e_h, ~, exit_flags(j)] = fsolve(g, 2, opts_fsolve);
    
    th_inf = acos(-1 / e_h);
    th_h = th_inf - alpha;
    rp_h = a_h * (1 - e_h^2);
    
    % Check Vincoli Fisici
    if th_inf > pi/2 && th_inf < pi && rp_h > R_earth
        
        h_ver = cross(rh_ver, r_inf) / norm(cross(rh_ver, r_inf)); 
        h_H_vec = sqrt(mu * a_h * (1 - e_h^2)) * h_ver; 
        
        control_orientamento = true;
        while control_orientamento
            
            N_H_vec = cross([0; 0; 1], h_H_vec); 
            N_H_vec = N_H_vec / norm(N_H_vec);
            i_H = acos(h_H_vec(3) / norm(h_H_vec));
            
            if N_H_vec(2) >= 0
                OM_H = acos(N_H_vec(1));
            else
                OM_H = 2*pi - acos(N_H_vec(1));
            end
            
            cos_beta = dot(N_H_vec, rh_ver);
            sin_beta = dot(cross(N_H_vec, rh_ver), (h_H_vec / norm(h_H_vec)));
            beta = atan2(sin_beta, cos_beta);
            om_H = beta - th_h;
            
            [~, v_hyp] = par2car(a_h, e_h, i_H, OM_H, om_H, th_h, mu);
            
            if dot(v_hyp, v_inf_earth_eci') < 0 
                h_H_vec = -h_H_vec; 

            else

                DV_values(j) = norm(v_hyp - v_park);
                e_h_vals(j)  = e_h;
                i_H_vals(j)  = i_H;
                OM_H_vals(j) = OM_H;
                om_H_vals(j) = om_H;
                
                control_orientamento = false;
            end
        end
    end
end


%% 4. ESTRAZIONE DELLA SOLUZIONE OTTIMALE (E DEL PERICENTRO)

% 4a. Soluzione Ottima (Minimo globale del Delta-V)
[min_DV, min_idx] = min(DV_values);
ottimo_th_i = theta_val(min_idx);

e_h_ott  = e_h_vals(min_idx);
i_H_ott  = i_H_vals(min_idx);
OM_H_ott = OM_H_vals(min_idx);
om_H_ott = om_H_vals(min_idx);

% Calcolo angoli geometrici per il report (Soluzione ottima)
[rh_vec, ~] = par2car(a_i, e_i, i_i, OM_i, om_i, ottimo_th_i, mu);
alpha_ott = acos(dot(rh_vec/norm(rh_vec), r_inf));
th_inf_ott = acos(-1 / e_h_ott);
th_h_ott = th_inf_ott - alpha_ott;
r_h_norm = norm(rh_vec);

deltaT = TOF_M_NoPrint(a_i, e_i, theta_iniziale, ottimo_th_i, mu);


% 4b. Soluzione al Pericentro (Indice 1, perché theta_val(1) = 0)
DV_pericentro = DV_values(1);
e_h_peri  = e_h_vals(1);
i_H_peri  = i_H_vals(1);
OM_H_peri = OM_H_vals(1);
om_H_peri = om_H_vals(1);

% Calcolo angoli geometrici per la soluzione al pericentro
[rh_vec, ~] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);
alpha_peri = acos(dot(rh_vec/norm(rh_vec), r_inf));
th_inf_peri = acos(-1 / e_h_peri);
th_h_peri = th_inf_peri - alpha_peri;


%% 5. STAMPA DEL REPORT DELLA MISSIONE
fprintf('\n=========================================================================\n');
fprintf('          MISSION ANALYSIS REPORT: EARTH DEPARTURE (SCENARIO 3)          \n');
fprintf('                      Non-Coplanar Escape Hyperbola                      \n');
fprintf('=========================================================================\n\n');

fprintf(' INIEZIONE AL PERICENTRO (Punto di Manovra)\n');
fprintf(' -------------------------------------------------------------------------\n');
fprintf('  Anomalia iniezione in orbita (th_in)   | %10.2f deg\n', rad2deg(0));
fprintf('  Semiasse maggiore (a_H)                | %10.2f km\n', a_h);
fprintf('  Eccentricita (e_H)                     | %10.6f\n', e_h_peri);
fprintf('  Inclinazione (i_H)                     | %10.2f deg\n', rad2deg(i_H_peri));
fprintf('  RAAN (OM_H)                            | %10.2f deg\n', rad2deg(OM_H_peri));
fprintf('  Arg. del Pericentro (om_H)             | %10.2f deg\n', rad2deg(om_H_peri+2*pi));
fprintf('  Anomalia asintotica (th_inf)           | %10.2f deg\n', rad2deg(th_inf_peri));
fprintf('  Anomalia sull''iperbole (th_H)          | %10.2f deg\n', rad2deg(th_h_peri));
fprintf('  Eccesso Iperbolico (V_inf)             | %10.4f km/s\n\n', v_inf_norm);

fprintf(' GEOMETRIA DELL''INIEZIONE (Punto di Manovra)\n');
fprintf(' -------------------------------------------------------------------------\n');
fprintf('  Anomalia iniezione in orbita (th_in)   | %10.2f deg\n', rad2deg(ottimo_th_i));
fprintf('  Distanza radiale dal centro (r_H)      | %10.2f km\n', r_h_norm);
fprintf('  Angolo di raccordo sull''iperbole (th_H)| %10.2f deg\n', rad2deg(th_h_ott));
fprintf('\n');

fprintf(' PARAMETRI IPERBOLE DI FUGA OTTIMIZZATA (Post-Manovra)\n');
fprintf(' -------------------------------------------------------------------------\n');
fprintf('  Semiasse maggiore (a_H)                | %10.2f km\n', a_h);
fprintf('  Eccentricita (e_H)                     | %10.6f\n', e_h_ott);
fprintf('  Inclinazione (i_H)                     | %10.2f deg\n', rad2deg(i_H_ott));
fprintf('  RAAN (OM_H)                            | %10.2f deg\n', rad2deg(OM_H_ott));
fprintf('  Arg. del Pericentro (om_H)             | %10.2f deg\n', rad2deg(om_H_ott+2*pi));
fprintf('  Anomalia asintotica (th_inf)           | %10.2f deg\n', rad2deg(th_inf_ott));
fprintf('  Eccesso Iperbolico (V_inf)             | %10.4f km/s\n\n', v_inf_norm);
fprintf('  Tempo di arrivo al punto di iniezione  | %10.2f s\n', deltaT);

fprintf(' COSTO DI MANOVRA (Delta-V)\n');
fprintf(' -------------------------------------------------------------------------\n');
fprintf('  Costo iniezione al pericentro          | %8.4f km/s\n', DV_pericentro);
fprintf('  Costo iniezione ottimizzata            | %8.4f km/s\n', min_DV);
fprintf('  ------------------------------------------------------------------------\n');
fprintf('  RISPARMIO NETTO OTTENUTO               | %8.4f km/s\n', abs(DV_pericentro - min_DV));

%% 6. PREPARAZIONE DATI PER IL PLOT 3D NON COMPLANARE
% Struttura per l'orbita di parcheggio
Park.a = a_i;
Park.e = e_i;
Park.i = i_i;
Park.OM = OM_i;
Park.om = om_i;

% Struttura per l'iperbole generata al Pericentro
Hyp_peri.a = a_h;
Hyp_peri.e = e_h_peri;       
Hyp_peri.i = i_H_peri;       
Hyp_peri.OM = OM_H_peri;     
Hyp_peri.om = om_H_peri;     
Hyp_peri.th_in = 0;          

% Struttura per l'iperbole ottima
Hyp_opt.a = a_h;
Hyp_opt.e = e_h_ott;
Hyp_opt.i = i_H_ott;         
Hyp_opt.OM = OM_H_ott;       
Hyp_opt.om = om_H_ott;       
Hyp_opt.th_in = ottimo_th_i;

scenery3_plot_noncoplanar(Park, Hyp_peri, Hyp_opt, mu, R_earth);