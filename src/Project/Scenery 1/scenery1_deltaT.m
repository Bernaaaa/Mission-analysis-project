% Scenery 1 - Delta T
% Transfer strategy minimizing the time of flight
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


DeltaT1 = TOF_M_NoPrint(a, e, th, pi, mu);
[DeltaV1BT, DeltaV2BT, tof_bit] = bitangentTransfer(a, e, a_f, e_f, 'ap', mu);

[dv_plane, om_post_plane, th_plane_change] = changeOrbitalPlaneDeltaT(a_f, e_f, i, OM, om, th, i_f, OM_f, mu);
DeltaT2 = TOF_M_NoPrint(a_f, e_f, 0, th_plane_change, mu);

[DeltaVCP, thi_omega_blt, thf_omega_blt] = changePericenterArg(a_f, e_f, om_post_plane, om_f, mu);
DeltaT3 = TOF_M_NoPrint(a_f, e_f, th_plane_change, thi_omega_blt(2), mu);
DeltaT4 = TOF_M_NoPrint(a_f, e_f, thf_omega_blt(2), th_f, mu);

DeltaT_TOT = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + tof_bit;
DeltaV_total = abs(DeltaV1BT) + abs(DeltaV2BT) + abs(dv_plane) + abs(DeltaVCP);

% flight time in days and hours and minutes and seconds
days = floor(DeltaT_TOT / (24*3600));
hours = floor(mod(DeltaT_TOT, 24*3600) / 3600);
minutes = floor(mod(DeltaT_TOT, 3600) / 60);
seconds = floor(mod(DeltaT_TOT, 60));


fprintf('\n=========================================================================\n');
fprintf('                 STRATEGY 2: MIN TOF MISSION SUMMARY                  \n');
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
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Bitangent Departure (BT-AP)', pi, abs(DeltaV1BT) + abs(DeltaV2BT));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Plane Change', th_plane_change, abs(dv_plane));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Arg. Pericenter Adj.', thi_omega_blt(2), abs(DeltaVCP));

% --- 3. SUMMARY FINALE ---
fprintf('\n=========================================================================\n');
fprintf('TOTALE DELTA-V          : %10.4f km/s\n', DeltaV_total);
fprintf('TEMPO TOTALE (TOF)      : %d d, %02d h, %02d m, %02d s\n', days, hours, minutes, seconds);
fprintf('=========================================================================\n\n');


% =========================================================================
% PREPARAZIONE DATI PER IL PLOT STANDARD (SCENARIO 1)
% =========================================================================

% 1. Initial Orbit (P1)
P1.a = a; P1.e = e; P1.i = i; P1.OM = OM; P1.om = om;
P1.th_in = th; P1.th_out = pi; % Arriva al apocentro
P1.name = 'Initial Orbit'; P1.maneuver = 'Start Position';

% 2. Bitangent Transfer (P2) - DA APOCENTRO A PERICENTRO
ra_initial = a*(1+e);
rp_final = a_f*(1-e_f);
P2.a = (ra_initial + rp_final)/2;       % Semi-major axis della bitangente di trasferimento 
P2.e = (ra_initial - rp_final)/(ra_initial + rp_final); % Eccentricità della bitangente di trasferimento
P2.i = i; P2.OM = OM; P2.om = om;
P2.th_in = pi; 
P2.th_out = 0; % Arriva al pericentro dell'orbita di trasferimento (th = 0)
P2.name = 'Bitangent Transfer (AP)'; P2.maneuver = '$\Delta V_{1BT}$';

% 3. Coasting post-Bitangente (P3)
P3.a = a_f; P3.e = e_f; P3.i = i; P3.OM = OM; P3.om = om;
P3.th_in = 0; 
P3.th_out = th_plane_change;
P3.name = 'Post-Bitangent Coasting'; P3.maneuver = '$\Delta V_{2BT}$ (Arrival)';

% 4. Coasting post-Cambio Piano (P4)
P4.a = a_f; P4.e = e_f; P4.i = i_f; P4.OM = OM_f; P4.om = om_post_plane;
P4.th_in = th_plane_change; 
P4.th_out = thi_omega_blt(2);
P4.name = 'Post-Plane Change Coasting'; P4.maneuver = '$\Delta V_{plane}$';

% 5. Final Orbit (P5)
P5.a = a_f; P5.e = e_f; P5.i = i_f; P5.OM = OM_f; P5.om = om_f;
P5.th_in = thf_omega_blt(2); 
P5.th_out = th_f;
P5.name = 'Final Orbit'; P5.maneuver = '$\Delta V_{\omega}$';

% --- Generazione Grafico ---
scenery1_plot_standard(P1, P2, P3, P4, P5);