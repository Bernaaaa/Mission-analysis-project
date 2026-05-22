close all; clear; clc;

%% --- SETUP PERCORSI UNIVERSALE ---
currentDir = fileparts(mfilename('fullpath'));

funcFolder = fullfile(currentDir, '..', '..', 'lab');

if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
else
    warning('Attenzione: La cartella delle funzioni non è stata trovata in: %s', funcFolder);
end

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


[DeltaVP_BLT, omf, theta] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);
Delta_t1 = TOF_M_NoPrint(a, e, th, theta, mu);

[DeltaVCP, thi, thf] = changePericenterArg(a, e, omf, om_f, mu);
Delta_t2 = TOF_M_NoPrint(a, e, theta, thi(2), mu);
Delta_t3 = TOF_M_NoPrint(a, e, thf(2), 0, mu);

[DeltaV1BT, DeltaV2BT, Delta_t4] = bitangentTransfer(a, e, a_f, e_f,'pa',mu);
Delta_t5 = TOF_M_NoPrint(a_f, e_f, pi, th_f, mu);

DeltaV_total = abs(DeltaVP_BLT) + abs(DeltaVCP) + abs(DeltaV1BT) + abs(DeltaV2BT);
DeltaT_TOT = abs(Delta_t1) + abs(Delta_t2) + abs(Delta_t3) + abs(Delta_t4) + abs(Delta_t5);


% flight time in days and hours and minutes and seconds
days = floor(DeltaT_TOT / (24*3600));
hours = floor(mod(DeltaT_TOT, 24*3600) / 3600);
minutes = floor(mod(DeltaT_TOT, 3600) / 60);
seconds = floor(mod(DeltaT_TOT, 60));


fprintf('\n=========================================================================\n');
fprintf('                 STRATEGY 1: STANDARD MISSION SUMMARY                  \n');
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
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Plane Change (Optimal DV)', theta, abs(DeltaVP_BLT));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Arg. Pericenter Adj.', thi(2), abs(DeltaVCP));
fprintf('%-30s | %-15.4f | %-12.4f\n', 'Bitangent Arrival (BT-PA)', pi, abs(DeltaV1BT) + abs(DeltaV2BT));

% --- 3. SUMMARY FINALE ---
fprintf('\n=========================================================================\n');
fprintf('TOTALE DELTA-V          : %10.4f km/s\n', DeltaV_total);
fprintf('TEMPO TOTALE (TOF)      : %d d, %02d h, %02d m, %02d s\n', days, hours, minutes, seconds);
fprintf('=========================================================================\n\n');

% =========================================================================
% PREPARAZIONE DATI PER IL PLOT STANDARD (SCENARIO 1)
% =========================================================================

% 1. Coasting Iniziale
P1.a = a; P1.e = e; P1.i = i; P1.OM = OM; P1.om = om;
P1.th_in = th; P1.th_out = theta;
P1.name = 'Initial Coasting'; P1.maneuver = 'Initial Pos';

% 2. Trans Orbit (Dopo Delta V Plane)
P2.a = a; P2.e = e; P2.i = i_f; P2.OM = OM_f; P2.om = omf;
P2.th_in = theta; P2.th_out = thi(2);
P2.name = 'Trans Orbit (Pre-CP)'; P2.maneuver = '$\Delta V_{plane}$';

% 3. Trans Orbit (Dopo Delta V CP)
P3.a = a; P3.e = e; P3.i = i_f; P3.OM = OM_f; P3.om = om_f;
P3.th_in = thf(2); P3.th_out = 2*pi; 
P3.name = 'Trans Orbit (Post-CP)'; P3.maneuver = '$\Delta V_{CP}$';

% 4. Trasferimento Bitangente
rp_current = a*(1-e); ra_f = a_f*(1+e_f);
P4.a = (rp_current + ra_f)/2; P4.e = (ra_f - rp_current)/(ra_f + rp_current);
P4.i = i_f; P4.OM = OM_f; P4.om = om_f;
P4.th_in = 0; P4.th_out = pi;
P4.name = 'Bitangent Transfer Leg'; P4.maneuver = '$\Delta V_{1BT}$';

% 5. Coasting su Orbita Finale
P5.a = a_f; P5.e = e_f; P5.i = i_f; P5.OM = OM_f; P5.om = om_f;
P5.th_in = pi; P5.th_out = th_f;
P5.name = 'Final Orbit Coasting'; P5.maneuver = '$\Delta V_{2BT}$';

% Lancia il plot
scenery1_plot_standard(P1, P2, P3, P4, P5);