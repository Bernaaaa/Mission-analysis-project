clear; clc; close all;

% --- 1. DATI DI INPUT
mu = 398600; a = 7000; e = 0.1; omi = deg2rad(30);
i_i = deg2rad(40); OMi = deg2rad(35);
i_f = deg2rad(20);  OMf = deg2rad(130.0);

% --- 2. ESECUZIONE DELLA FUNZIONE ---
[DeltaV, omf, theta] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu);

% --- 3. NORMALIZZAZIONE PER IL TEST ---
theta_deg = mod(rad2deg(theta), 360);
omf_deg   = mod(rad2deg(omf), 360);
theta_rad = deg2rad(theta_deg);
omf_rad   = deg2rad(omf_deg);

% --- 4. VERIFICA MATEMATICA ---
getPos = @(i, OM, om, th) (a*(1-e^2)/(1+e*cos(th))) * [ ...
    cos(OM)*cos(om+th) - sin(OM)*sin(om+th)*cos(i); ...
    sin(OM)*cos(om+th) + cos(OM)*sin(om+th)*cos(i); ...
    sin(om+th)*sin(i) ];

P_inizio = getPos(i_i, OMi, omi, theta_rad);
P_fine   = getPos(i_f, OMf, omf_rad, theta_rad);
distanza_errore = norm(P_inizio - P_fine);

% --- 5. PRINT RISULTATI ---
fprintf('--- REPORT VALIDAZIONE ---\n');
fprintf('Punto di Manovra (Theta): %.2f deg\n', theta_deg);
fprintf('Nuovo Pericentro (omf):   %.2f deg\n', omf_deg);
fprintf('Costo Manovra (DeltaV):   %.4f km/s\n', DeltaV);
fprintf('Errore di Coincidenza:    %.6e km\n', distanza_errore);
if distanza_errore < 1e-3
    fprintf('RISULTATO: PERFETTO. Il punto è un''intersezione reale.\n');
else
    fprintf('RISULTATO: ERRORE. Il punto non tocca entrambe le orbite.\n');
end
fprintf('--------------------------\n');

% --- 6. PLOT 3D ---
th_vec = linspace(0, 2*pi, 500);
orb_i = zeros(3, 500); orb_f = zeros(3, 500);
for j = 1:500
    orb_i(:,j) = getPos(i_i, OMi, omi, th_vec(j));
    orb_f(:,j) = getPos(i_f, OMf, omf_rad, th_vec(j));
end

figure('Color', 'w', 'Name', 'Validazione Orbitale');
plot3(orb_i(1,:), orb_i(2,:), orb_i(3,:), 'b--', 'LineWidth', 1.2); hold on;
plot3(orb_f(1,:), orb_f(2,:), orb_f(3,:), 'r-', 'LineWidth', 1.5);
plot3(P_inizio(1), P_inizio(2), P_inizio(3), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 10);

% Abbellimenti
grid on; axis equal; xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend('Orbita Iniziale', 'Orbita Finale', 'Punto di Manovra');
title('Verifica Geometrica Cambio di Piano');
view(135, 30);