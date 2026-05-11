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


%% STRATEGY 2: Bielliptic Transfer + Plane Change at Apogee
fprintf('\n--- STRATEGY 2: Bielliptic Transfer + Plane Change at Apogee ---\n');

% --- SETTAGGIO GRIGLIA DI RICERCA ---
ra_range = linspace(100000, 350000, 500); % 50 valori di raggio apocentro
e_range = linspace(0.01, 0.95, 500);      % 50 valori di eccentricità

best_dv = inf;
best_ra = 0;
best_e = 0;

fprintf('Ricerca della coppia ottimale in corso...\n');

for r_test = ra_range
    for e_test = e_range
        
        % 1. Calcolo semiasse temporaneo
        a_test = r_test / (1 + e_test);
        
        % 2. Salta se esce dalla SOI o se l'orbita non è fisica
        rp_test = a_test * (1 - e_test);
        if rp_test < 6378 
            continue;
        end
        
        try 
            [DV1a, DV1b] = bitangentTransfer(a, e, a_test, e_test, 'pa', mu);
            [DVP, om_p] = changeOrbitalPlaneSpecial(a_test, e_test, i, OM, om, i_f, OM_f, mu);
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

fprintf('TROVATO! Miglior ra_t: %.2f km, Miglior e_t: %.4f\n\n', best_ra, best_e);

ra_t = best_ra; % Apocenter of the transfer orbit (chosen to be very high for a more efficient bielliptic transfer)
e_t = best_e; % Eccentricity of the transfer orbit
a_t = ra_t / (1 + e_t); % Semi-major axis of the transfer orbit

[DeltaV1BT, DeltaV2BT, Delta_t_bt] = bitangentTransfer(a, e, a_t, e_t, 'aa', mu);
Delta_T1 = TOF_M(a, e, th, 0, mu); % Time from initial orbit to apogee of transfer orbit

fprintf('DeltaV for first leg of bielliptic transfer: %.4f km/s\n', abs(DeltaV1BT) + abs(DeltaV2BT));

[DeltaVP_BLT, om_post_plane_blt, theta_plane_blt] = changeOrbitalPlane(a_t, e_t, i, OM, om, i_f, OM_f, mu); % Plane change at apogee of transfer orbit
Delta_T2 = TOF_M(a_t, e_t, pi, theta_plane_blt, mu); % Time from apogee of transfer orbit to plane change point


[DeltaV_omega_blt, thi_omega_blt, thf_omega_blt] = changePericenterArg(a_t, e_t, om_post_plane_blt, om_f, mu); % Clean-up at final orbit
Delta_T3 = TOF_M(a_t, e_t, theta_plane_blt, thi_omega_blt(2), mu); % Time from plane change point to pericenter argument change point
Delta_T4 = TOF_M(a_t, e_t, thf_omega_blt(2), 0, mu); % Time from pericenter argument change point to final orbit

fprintf('DeltaV for argument of pericenter change: %.4f km/s\n', abs(DeltaV_omega_blt));

[DeltaV3BT, DeltaV4BT, Delta_t_bt2] = bitangentTransfer(a_t, e_t, a_f, e_f, 'ap', mu); % Final transfer to target orbit
Delta_T5 = TOF_M(a_t, e_t, 0, pi, mu); % Time from pericenter of transfer orbit to apogee (where the second leg of the bielliptic transfer happens)

fprintf('DeltaV for second leg of bielliptic transfer: %.4f km/s\n', abs(DeltaV3BT) + abs(DeltaV4BT));

DeltaV_Totale_BLT = abs(DeltaV1BT) + abs(DeltaV2BT) + abs(DeltaVP_BLT) + abs(DeltaV_omega_blt) + abs(DeltaV3BT) + abs(DeltaV4BT);
fprintf('Total Bielliptic Transfer DeltaV: %.4f km/s\n', DeltaV_Totale_BLT);

DeltaT_TOT = Delta_T1 + Delta_T2 + Delta_T3 + Delta_T4 + Delta_T5 + Delta_t_bt + Delta_t_bt2; % Total time of flight for the bielliptic transfer
%flight time in DD HH MM SS format
days = floor(DeltaT_TOT / (24*3600));
hours = floor(mod(DeltaT_TOT, 24*3600) / 3600);
minutes = floor(mod(DeltaT_TOT, 3600) / 60);
seconds = floor(mod(DeltaT_TOT, 60));

fprintf('Total Time of Flight for Bielliptic Transfer: %d days, %d hours, %d minutes, %d seconds\n', days, hours, minutes, seconds);
