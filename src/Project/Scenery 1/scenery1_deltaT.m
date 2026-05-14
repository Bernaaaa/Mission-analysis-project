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

fprintf('Initial parameters\n a = %.2f\n e = %.6f\n i = %.6f\n OM = %.6f\n om = %.6f\n th = %.6f\n\n', a, e, i, OM, om, rad2deg(th));
fprintf('Final parameters\n af = %.2f\n ef = %.6f\n if = %.6f\n OMf = %.6f\n omf = %.6f\n thf = %.6f\n\n', a_f, e_f, i_f, OM_f, om_f, rad2deg(th_f));

DeltaT1 = TOF_M_NoPrint(a, e, th, 0, mu);
[DV1, DV2, tof_bit] = bitangentTransfer(a, e, a_f, e_f, 'pp', mu);


[dv_plane, om_post_plane, th_plane_change] = changeOrbitalPlaneDeltaT(a_f, e_f, i, OM, om, th, i_f, OM_f, mu);
DeltaT2 = TOF_M_NoPrint(a_f, e_f, 0, th_plane_change, mu);


[DV3, thi_omega_blt, thf_omega_blt] = changePericenterArg(a_f, e_f, om_post_plane, om_f, mu);
DeltaT3 = TOF_M_NoPrint(a_f, e_f, th_plane_change, thi_omega_blt(2), mu);
DeltaT4 = TOF_M_NoPrint(a_f, e_f, thf_omega_blt(2), th_f, mu);

DeltaT_TOT = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + tof_bit;

fprintf('\nTotal maneuver cost deltaV = %.4f km/s\n\n', abs(DV1) + abs(DV2) + abs(dv_plane) + abs(DV3));

% flight time in days and hours and minutes and seconds
days = floor(DeltaT_TOT / (24*3600));
hours = floor(mod(DeltaT_TOT, 24*3600) / 3600);
minutes = floor(mod(DeltaT_TOT, 3600) / 60);
seconds = floor(mod(DeltaT_TOT, 60));

fprintf('Total Time of Flight for transfer: %d days, %d hours, %d minutes, %d seconds\n', days, hours, minutes, seconds);