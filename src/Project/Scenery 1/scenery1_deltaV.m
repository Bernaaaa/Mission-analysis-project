%Strategy 2: Bielliptic Transfer + Plane Change at Apogee + Clean-up at Final Orbit
% The objective is to maximize the efficiency of the transfer by performing the plane change at the apogee of a 
% highly elliptical transfer orbit, where the velocity is lower, and then doing a final clean-up maneuver 
% to adjust the argument of pericenter.
close all; clear; clc;

%% --- SETUP PERCORSI UNIVERSALE ---
currentDir = fileparts(mfilename('fullpath'));

funcFolder = fullfile(currentDir, '..', '..', 'lab');

if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
else
    warning('Attenzione: La cartella delle funzioni non è stata trovata in: %s', funcFolder);
end

%% --- PARAMETRI INIZIALI ---

mu = 398600.4418;

a = 24400.00;
e = 0.728300;
i = 0.104700;
OM = 2.361000;
om = 3.107000;
th = 2.135000;

r_xf = -7090.590200;
r_yf = -5612.557300;
r_zf = 3948.902900;
v_xf = 5.698000;
v_yf = -5.995000;
v_zf = 1.710000;

rf = [r_xf; r_yf; r_zf];
vf = [v_xf; v_yf; v_zf];

rfm = norm(rf,2);
vfm = norm(vf,2);

[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rf, vf, mu);

ra_range = linspace(50000, 315000, 500);
e_range = linspace(0.01, 0.95, 500);      

best_dv = inf;
best_ra = 0;
best_e = 0;

fprintf('Ricerca della coppia ottimale in corso...  ');

for r_test = ra_range
    for e_test = e_range
        
        a_test = r_test / (1 + e_test);
        
        rp_test = a_test * (1 - e_test);
        if rp_test < 6578 
            continue;
        end
        
        try 
            [DV1a, DV1b] = bitangentTransfer(a, e, a_test, e_test, 'pa', mu);
            [DVP, om_p] = changeOrbitalPlaneNoPrint(a_test, e_test, i, OM, om, i_f, OM_f, mu);
            [DVom] = changePericenterArg(a_test, e_test, om_p, om_f, mu);
            [DV2a, DV2b] = bitangentTransfer(a_test, e_test, a_f, e_f, 'ap', mu);
            
            DV_tot = abs(DV1a) + abs(DV1b) + abs(DVP) + abs(DVom) + abs(DV2a) + abs(DV2b);
            
            % --- Salvataggio del Record ---
            if DV_tot < best_dv
                best_dv = DV_tot;
                best_ra = r_test;
                best_e = e_test;
            end
        catch
            continue; % Se la geometria è impossibile, passa alla prossima
        end
    end
end

fprintf('TROVATO!\n');

ra_t = best_ra; % Apocenter of the transfer orbit (chosen to be very high for a more efficient bielliptic transfer)
e_t = best_e; % Eccentricity of the transfer orbit
a_t = ra_t / (1 + e_t); % Semi-major axis of the transfer orbit

[DeltaV1BT, DeltaV2BT, Delta_t_bt] = bitangentTransfer(a, e, a_t, e_t, 'pa', mu);
Delta_T1 = TOF_M_NoPrint(a, e, th, 0, mu);                                                                          % Time from initial orbit to apogee of transfer orbit

[DeltaVP_BLT, om_post_plane_blt, theta_plane_blt] = changeOrbitalPlane(a_t, e_t, i, OM, om, i_f, OM_f, mu);         % Plane change at apogee of transfer orbit
Delta_T2 = TOF_M_NoPrint(a_t, e_t, pi, theta_plane_blt, mu);                                                        % Time from apogee of transfer orbit to plane change point

[DeltaV_omega_blt, thi_omega_blt, thf_omega_blt] = changePericenterArg(a_t, e_t, om_post_plane_blt, om_f, mu);      % Clean-up at final orbit
Delta_T3 = TOF_M_NoPrint(a_t, e_t, theta_plane_blt, thi_omega_blt(1), mu);                                          % Time from plane change point to pericenter argument change point
Delta_T4 = TOF_M_NoPrint(a_t, e_t, thf_omega_blt(1), pi, mu);                                                       % Time from pericenter argument change point to final orbit

[DeltaV3BT, DeltaV4BT, Delta_t_bt2] = bitangentTransfer(a_t, e_t, a_f, e_f, 'ap', mu);                              % Final transfer to target orbit
Delta_T5 = TOF_M_NoPrint(a_f, e_f, 0, th_f, mu);                                                                    % Time from pericenter argument change point to final orbit


DeltaV_Totale_BLT = abs(DeltaV1BT) + abs(DeltaV2BT) + abs(DeltaVP_BLT) + abs(DeltaV_omega_blt) + abs(DeltaV3BT) + abs(DeltaV4BT);
DeltaT_TOT = Delta_T1 + Delta_T2 + Delta_T3 + Delta_T4 + Delta_T5 + Delta_t_bt + Delta_t_bt2;                       % Total time of flight for the bielliptic transfer

%flight time in DD HH MM SS format
days = floor(DeltaT_TOT / (24*3600));
hours = floor(mod(DeltaT_TOT, 24*3600) / 3600);
minutes = floor(mod(DeltaT_TOT, 3600) / 60);
seconds = floor(mod(DeltaT_TOT, 60));

fprintf('\n=========================================================================\n');
fprintf('                 STRATEGY 2: BIELLIPTIC MISSION SUMMARY                  \n');
fprintf('=========================================================================\n');

% --- 1. COMPARAZIONE ORBITALE ---
fprintf('\n%-20s | %-15s | %-15s\n', 'Parametro', 'Orbite Iniziale', 'Orbite Finale');
fprintf('-------------------------------------------------------------------------\n');
fprintf('%-20s | %-15.2f | %-15.2f\n', 'Semi-asse (a)', a, a_f);
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Eccentricità (e)', e, e_f);
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Inclinazione (i)', i, i_f);
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Nodo (OM)', OM, OM_f);
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Pericentro (om)', om, om_f);
fprintf('%-20s | %-15.4f | %-15.4f\n', 'Anomalia (th)', th, th_f);

% --- 2. LOG MANOVRE ---
fprintf('\n%-30s | %-15s | %-12s\n', 'Manovra', 'Anomalia (rad)', 'DeltaV (km/s)');
fprintf('-------------------------------------------------------------------------\n');
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Bitangent Departure (BT1-PA)', 0, abs(DeltaV1BT) + abs(DeltaV2BT));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Plane Change (Optimal DV)', theta_plane_blt, abs(DeltaVP_BLT));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Arg. Pericenter Adj.', thi_omega_blt(1), abs(DeltaV_omega_blt));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Bitangent Arrival (BT2-AP)', 0, abs(DeltaV3BT) + abs(DeltaV4BT));

% --- 3. SUMMARY FINALE ---
fprintf('\n=========================================================================\n');
fprintf('TOTALE DELTA-V          : %10.4f km/s\n', DeltaV_Totale_BLT);
fprintf('TEMPO TOTALE (TOF)      : %d d, %02d h, %02d m, %02d s\n', days, hours, minutes, seconds);
fprintf('=========================================================================\n\n');

% --- PREPARAZIONE DATI PER IL PLOT ---
GTO.a = a; GTO.e = e; GTO.i = i; GTO.OM = OM; GTO.om = om; GTO.th = th; GTO.mu = mu;
Park.a = a_f; Park.e = e_f; Park.i = i_f; Park.OM = OM_f; Park.om = om_f; Park.th = th_f; Park.mu = mu;

% Ramo 1: punti esatti delle manovre interne
Trans1.a = a_t; Trans1.e = e_t; Trans1.i = i; Trans1.OM = OM; Trans1.om = om; Trans1.mu = mu;
Trans1.theta_plane = theta_plane_blt; Trans1.om_mid = om_post_plane_blt;
Trans1.thi_omega = thi_omega_blt(1);

% Ramo 2: punto di uscita dalla manovra di omega
Trans2.a = a_t; Trans2.e = e_t; Trans2.i = i_f; Trans2.OM = OM_f; Trans2.om = om_f; Trans2.mu = mu;
Trans2.thf_omega = thf_omega_blt(1);

scenery1_plot_deltaV(GTO, Park, Trans1, Trans2);