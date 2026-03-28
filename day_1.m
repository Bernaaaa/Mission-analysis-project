clc 
clear 
close all

%% LEZIONE 1 - 27/03/2026 - SCENARIO 1
%DATI LABORATORIO

%Orbita iniziale
 a_i = 24400.00;
 e_i = 0.728300;
 incl_i = 0.104700;
 OM_i = 2.361000;
 om_i = 3.107000;
 theta_i = 2.135000;

 %Orbita finale
  rx_f = -7090.590200;
  rv_f = -5612.557300;
  rz_f =3948.902900;
  vx_f = 5.698000;
  vy_f = -5.995000;
  vz_f = 1.710000;

  %Costanti gravitazionali
  mu_earth = 398600; %km^3/s^2
  mu_sun = 132712440017.99; %km^3/s^2


%completiamo i dati mancanti dell'orbita finale
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par([rx_f; rv_f; rz_f], [vx_f; vy_f; vz_f], mu_earth);

%completiamo i dati mancanti dell'orbita iniziale
[rr, vv] = par2car(a_i, e_i, incl_i, OM_i, om_i, theta_i, mu_earth);

rx_i = rr(1);
ry_i = rr(2);
rz_i = rr(3);

vx_i = vv(1);
vy_i = vv(2);
vz_i = vv(3);

%plot dell'orbita iniziale
figure('Name', 'Orbita iniziale')
plotOrbit(a_i, e_i, incl_i, OM_i, om_i, 0, 2*pi, 0.01, mu_earth)

%plot dell'orbita finale
figure('Name', 'Orbita finale')
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu_earth)
