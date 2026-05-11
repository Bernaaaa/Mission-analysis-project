%a_i [km] e_i inc_i [rad] RAAN_i [rad] w_i [rad] theta_i [rad] rx_f [km] ry_f [km] rz_f [km] vx_f [km/s] vy_f [km/s] vz_f [km/s]
clear; close all; clc;

mu = 398600.4418;
a = 24400.00;
e = 0.728300;
i = 0.104700;
OM = 1.289000;
om = 3.107000;
th = 0.510900;

r_xf = -6103.007500;
r_yf = 13604.661000;
r_zf = 7991.352600;

v_xf = -5.163000;
v_yf = -3.151000;
v_zf = 1.421000;


rf = [r_xf; r_yf; r_zf];
vf = [v_xf; v_yf; v_zf];

rfm = norm(rf,2);
vfm = norm(vf,2);


[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rf, vf, mu);

rad2deg(om_f)

fprintf('Final parameters\n af = %.2f\n ef = %.6f\n if = %.6f\n OMf = %.6f\n omf = %.6f\n thf = %.6f\n\n', a_f, e_f, i_f, OM_f, om_f, th_f);


[DeltaV1BT, DeltaV2BT, Delta_t4] = bitangentTransfer(a, e, a_f, e_f,'pa',mu);
fprintf('DeltaV for bitangent transfer: %.4f km/s\n', DeltaV1BT + DeltaV2BT);
Delta_t5 = TOF_M(a, e, th, 0, mu);


[DeltaV, omf, theta] = changeOrbitalPlane(a_f, e_f, i, OM, om, i_f, OM_f, mu);
Delta_t1 = TOF_M(a_f, e_f, th, theta, mu);

[DeltaVCP, thi, thf] = changePericenterArg(a_f, e_f, omf-pi, om_f, mu);
fprintf('DeltaV for argument of pericenter change: %.4f km/s\n', DeltaVCP);
Delta_t2 = TOF_M(a_f, e_f, theta, thi(1), mu);
Delta_t3 = TOF_M(a_f, e_f, thf(1), th_f, mu);



DeltaV_total = abs(DeltaV) + abs(DeltaVCP) + abs(DeltaV1BT) + abs(DeltaV2BT);
DeltaT_total = abs(Delta_t1) + abs(Delta_t2) + abs(Delta_t3) + abs(Delta_t4) + abs(Delta_t5);

fprintf('Total DeltaV: %.4f km/s\n', DeltaV_total);
fprintf('Total Time of Flight: %.2f seconds\n', DeltaT_total);