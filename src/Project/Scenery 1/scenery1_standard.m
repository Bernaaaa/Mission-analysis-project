close all; clear; clc;

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

fprintf('Final parameters\n af = %.2f\n ef = %.6f\n if = %.6f\n OMf = %.6f\n omf = %.6f\n thf = %.6f\n\n', a_f, e_f, i_f, OM_f, om_f, th_f);

[DeltaV, omf, theta] = changeOrbitalPlane(a, e, i, OM, om, i_f, OM_f, mu);
Delta_t1 = TOF_M(a, e, th, theta, mu);

[DeltaVCP, thi, thf] = changePericenterArg(a, e, omf, om_f, mu);
fprintf('DeltaV for argument of pericenter change: %.4f km/s\n', DeltaVCP);
Delta_t2 = TOF_M(a, e, theta, thi(2), mu);
Delta_t3 = TOF_M(a, e, thf(2), 0, mu);

[DeltaV1BT, DeltaV2BT, Delta_t4] = bitangentTransfer(a, e, a_f, e_f,'pa',mu);
fprintf('DeltaV for bitangent transfer: %.4f km/s\n', DeltaV1BT + DeltaV2BT);
Delta_t5 = TOF_M(a_f, e_f, pi, th_f, mu);

DeltaV_total = abs(DeltaV) + abs(DeltaVCP) + abs(DeltaV1BT) + abs(DeltaV2BT);
DeltaT_TOT = abs(Delta_t1) + abs(Delta_t2) + abs(Delta_t3) + abs(Delta_t4) + abs(Delta_t5);

fprintf('Total DeltaV: %.4f km/s\n', DeltaV_total);

% flight time in days and hours and minutes and seconds
days = floor(DeltaT_TOT / (24*3600));
hours = floor(mod(DeltaT_TOT, 24*3600) / 3600);
minutes = floor(mod(DeltaT_TOT, 3600) / 60);
seconds = floor(mod(DeltaT_TOT, 60));

fprintf('Total Time of Flight for transfer: %d days, %d hours, %d minutes, %d seconds\n', days, hours, minutes, seconds);

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