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

% trying more economic transfers

% 1. CAMBIO FORMA (Bitangent Transfer 'pa')
% Dal pericentro dell'orbita iniziale all'apocentro dell'orbita finale.
[DeltaV1BT, DeltaV2BT, Delta_t_bt] = bitangentTransfer(a, e, a_f, e_f, 'pa', mu);
fprintf('DeltaV for bitangent transfer: %.4f km/s\n', abs(DeltaV1BT) + abs(DeltaV2BT));

% 2. CAMBIO PIANO 
% Ora sei sull'orbita FINALE (a_f, e_f). Esegui il cambio di piano qui.
% Nota: passiamo a_f ed e_f, non numeri inventati!
[DeltaVP, om_post_plane, theta_plane] = changeOrbitalPlane(a_f, e_f, i, OM, om, i_f, OM_f, mu);
fprintf('DeltaV for plane change: %.4f km/s\n', abs(DeltaVP));

% 3. CAMBIO PERICENTRO (Clean-up)
% Sei sempre sull'orbita finale (a_f, e_f), ma l'argomento è diventato om_post_plane.
[DeltaV_omega, thi_omega, thf_omega] = changePericenterArg(a_f, e_f, om_post_plane, om_f, mu);
fprintf('DeltaV for argument of pericenter change: %.4f km/s\n', abs(DeltaV_omega));

% --- TOTALI ---
DeltaV_Totale_Ottimizzato = abs(DeltaV1BT) + abs(DeltaV2BT) + abs(DeltaVP) + abs(DeltaV_omega);
fprintf('\nTotal Optimized DeltaV: %.4f km/s\n', DeltaV_Totale_Ottimizzato);

