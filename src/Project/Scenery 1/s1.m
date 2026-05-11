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
DeltaT_total = abs(Delta_t1) + abs(Delta_t2) + abs(Delta_t3) + abs(Delta_t4) + abs(Delta_t5);

fprintf('Total DeltaV: %.4f km/s\n', DeltaV_total);
fprintf('Total Time of Flight: %.2f seconds\n', DeltaT_total);