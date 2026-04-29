
% --- 1. DATI DI INPUT
mu = 398600; a = 12000; e = 0.2; omi = deg2rad(310);
i_i = deg2rad(45); OMi = deg2rad(80);
i_f = deg2rad(25);  OMf = deg2rad(55.0);

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

% --- 6. RENDER 3D "SPACEX MISSION CONTROL - 16K" ---
res = 16000; % Risoluzione vettoriale per orbite (16K per massimizzare la fluidità e dettaglio)
th_vec = linspace(0, 2*pi, res);
orb_i = zeros(3, res); orb_f = zeros(3, res);
for j = 1:res
    orb_i(:,j) = getPos(i_i, OMi, omi, th_vec(j));
    orb_f(:,j) = getPos(i_f, OMf, omf_rad, th_vec(j));
end

if exist('fig', 'var') && ishandle(fig), clf(fig); else, fig = figure; end
set(fig, 'Color', [0.03 0.03 0.04], 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.85]);
ax = gca; hold on; axis equal; grid off;
set(ax, 'Color', [0.03 0.03 0.04], 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
view(145, 25);
set(fig, 'WindowState', 'maximized');




% --- 6a. TERRA 16K (High-Contrast Topography) ---
r_pericentro = a * (1 - e); 
R_earth = 6371;  
[X_e, Y_e, Z_e] = sphere(300); 
load topo;

topo_mask = topo > -500; 
topo_tex = zeros(size(topo,1), size(topo,2), 3);
ocean_c = [0.05 0.07 0.1];  
land_c  = [0.25 0.3 0.35];  

for i = 1:3
    layer = topo_tex(:,:,i);
    layer(~topo_mask) = ocean_c(i);
    layer(topo_mask) = land_c(i);
    topo_tex(:,:,i) = layer;
end
surface(X_e*R_earth, Y_e*R_earth, Z_e*R_earth, topo_tex, ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'AmbientStrength', 0.8);

surface(X_e*R_earth*1.02, Y_e*R_earth*1.02, Z_e*R_earth*1.02, ...
    'FaceColor', [0.2 0.5 1], 'EdgeColor', 'none', 'FaceAlpha', 0.05);

% --- 6b. ORBITE "PLASMA-GLOW" (16K NEON RENDERING) ---
c_white = [1 1 1];       % White Core
c_red   = [1 0.05 0.05]; % Crimson Red

% --- RENDERING ORBITA INIZIALE (NEON WHITE) ---
% Layer 1: Core continuo
h_orb_i = plot3(orb_i(1,:), orb_i(2,:), orb_i(3,:), 'Color', c_white, 'LineWidth', 1.4, ...
    'LineStyle', '-', 'DisplayName', 'INITIAL ORBIT STABLE');

% --- RENDERING ORBITA FINALE (PULSE RED) ---
% Layer 1: Core continuo 
h_orb_f = plot3(orb_f(1,:), orb_f(2,:), orb_f(3,:), 'Color', c_red, 'LineWidth', 1.6, ...
    'LineStyle', '-', 'DisplayName', 'TARGET PLANE LOCK'); % <--- Ora è una linea continua

    

% --- 6c. ECI VECTORS (Frecce e Testi Dinamici Adattivi) ---
r_apocentro = a * (1 + e);
axis_len = r_apocentro * 1.2; % Gli assi si allungano sempre il 20% oltre il punto più lontano dell'orbita

% Griglia polare di base (ora si adatta dinamicamente)
th_grid = linspace(0, 2*pi, 200);
for r_step = linspace(R_earth*1.2, axis_len*0.8, 4)
    plot3(r_step*cos(th_grid), r_step*sin(th_grid), zeros(size(th_grid)), ...
        'Color', [1 1 1 0.03], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Funzione per disegnare vettori con punte a freccia
draw_arrow = @(p1, p2, col, w) quiver3(p1(1), p1(2), p1(3), p2(1)-p1(1), p2(2)-p1(2), p2(3)-p1(3), ...
    0, 'Color', col, 'LineWidth', w, 'MaxHeadSize', 0.08);

% Disegno dei vettori con lunghezza dinamica
draw_arrow([0 0 0], [axis_len 0 0], [0.2 0.6 1], 1.5); % Vernal X
draw_arrow([0 0 0], [0 axis_len 0], [0.2 0.6 1 0.3], 1); % Y
draw_arrow([0 0 0], [0 0 axis_len], [0.8 0.8 0.8], 1.5); % Polar Z

% Posizionamento delle etichette (si spostano sempre un 5% oltre la punta della freccia)
text(axis_len*1.05, 0, 0, 'X_{VERNAL}', 'Color', [0.2 0.6 1], 'FontName', 'Helvetica', 'FontSize', 10);
text(0, 0, axis_len*1.05, 'Z_{POLAR}', 'Color', [0.8 0.8 0.8], 'FontName', 'Helvetica', 'FontSize', 10);
% Asse Y (Completa la terna ECI)
text(0, axis_len*1.05, 0, 'Y_{ECI}', 'Color', [0.2 0.6 1], 'FontName', 'Consolas', 'FontSize', 10);

% --- 6d. MANEUVER TARGETING ---
c_amber = [1.0, 0.8, 0.0];

% Vettore posizionale dal centro della Terra
line([0 P_inizio(1)], [0 P_inizio(2)], [0 P_inizio(3)], 'Color', [1 0.2 0 0.2], 'LineStyle', '-', 'LineWidth', 1);

% --- Colore Verde ---
c_green = [0.0, 0.8, 0.3]; 

% --- TARGET LOCK (Mirino Tactical Green) ---
% Layer 1: Punto centrale
plot3(P_inizio(1), P_inizio(2), P_inizio(3), 'go', 'MarkerSize', 4, 'MarkerFaceColor', 'g'); 

% Layer 3: Croce di puntamento
plot3(P_inizio(1), P_inizio(2), P_inizio(3), '+', 'Color', c_green, 'MarkerSize', 12, 'LineWidth', 0.8);

% Layer 4: Box angolare
plot3(P_inizio(1), P_inizio(2), P_inizio(3), 's', 'Color', c_green, 'MarkerSize', 15, 'LineWidth', 0.5);
% ==================== CALCOLO VETTORE DI SPINTA REALE ====================
p_val = a * (1 - e^2); h_val = sqrt(mu * p_val);
u_i = omi + theta_rad;
v_i = (mu/h_val) * [-cos(OMi)*(sin(u_i)+e*sin(omi))-sin(OMi)*(cos(u_i)+e*cos(omi))*cos(i_i); -sin(OMi)*(sin(u_i)+e*sin(omi))+cos(OMi)*(cos(u_i)+e*cos(omi))*cos(i_i); (cos(u_i)+e*cos(omi))*sin(i_i)];
u_f = omf_rad + theta_rad;
v_f = (mu/h_val) * [-cos(OMf)*(sin(u_f)+e*sin(omf_rad))-sin(OMf)*(cos(u_f)+e*cos(omf_rad))*cos(i_f); -sin(OMf)*(sin(u_f)+e*sin(omf_rad))+cos(OMf)*(cos(u_f)+e*cos(omf_rad))*cos(i_f); (cos(u_f)+e*cos(omf_rad))*sin(i_f)];

deltaV_vec = v_f - v_i;
arrow_len = 3000; 
thrust_v = P_inizio + (deltaV_vec / norm(deltaV_vec)) * arrow_len;
% =========================================================================

% DISEGNO VETTORE
quiver3(P_inizio(1), P_inizio(2), P_inizio(3), ...
    thrust_v(1)-P_inizio(1), thrust_v(2)-P_inizio(2), thrust_v(3)-P_inizio(3), ...
    0, 'Color', c_amber, 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale', 'off');

text(thrust_v(1), thrust_v(2), thrust_v(3)+800, ...
    'F_{THRUST}', ... 
    'Color', c_amber, 'FontName', 'Consolas', 'FontSize', 10, 'FontWeight', 'bold', ...
    'Interpreter', 'tex', 'HorizontalAlignment', 'right');


    % --- 6i. ORBITAL DIRECTION ANALYSIS (Prograde vs Retrograde) ---
% Calcolo momento angolare
h_vec = cross(P_inizio, v_i); 
is_prograde = h_vec(3) > 0; % Se la componente Z è positiva, è antioraria

if is_prograde
    rot_str = 'PROGRADE (CCW)';
    rot_col = [0.2, 0.8, 1.0]; % Ciano
    arc_dir = 1; % Direzione per il disegno dell'arco
else
    rot_str = 'RETROGRADE (CW)';
    rot_col = [1.0, 0.4, 0.2]; % Arancio/Rosso
    arc_dir = -1;
end

% --- DISEGNO INDICATORE DI ROTAZIONE (Ring Arrow attorno alla Terra) ---
r_ring = R_earth * 1.3; % Raggio dell'indicatore poco sopra la superficie
phi_ring = linspace(0, arc_dir * pi/2, 50); % Un quarto di cerchio per indicare il verso
x_ring = r_ring * cos(phi_ring);
y_ring = r_ring * sin(phi_ring);
z_ring = zeros(size(phi_ring));

plot3(x_ring, y_ring, z_ring, 'Color', [rot_col 0.5], 'LineWidth', 2, 'LineStyle', '--');

% Freccia sulla punta dell'arco per indicare il verso
quiver3(x_ring(end), y_ring(end), z_ring(end), ...
    -arc_dir * y_ring(end)*0.3, arc_dir * x_ring(end)*0.3, 0, ...
    0, 'Color', rot_col, 'LineWidth', 2, 'MaxHeadSize', 2);


% Stampa log
fprintf('Assetto Orbitale: %s\n', rot_str);

% 6e. TELEMETRIA DINAMICA - Posizionamento dinamico del testo
% Dati Telemetrici Estesi al Punto di Manovra
P_text = P_inizio * 2.5; % Lo allontaniamo un po' di più per dare respiro

% 1. Linea di richiamo con trasparenza e ordine di profondità gestito
line([P_inizio(1) P_text(1)], [P_inizio(2) P_text(2)], [P_inizio(3) P_text(3)], ...
    'Color', [0.7 0.7 0.7 0.3], 'LineStyle', ':', 'LineWidth', 1.2, ...
    'Clipping', 'on'); % Si nasconde se finisce dietro la Terra

% 2. Textbox con gestione avanzata dei layer
t_hud = text(P_text(1), P_text(2), P_text(3), ...
    {'\bf[ BURN COORD LOCK ]', ...
     sprintf('\\rm{THETA     : %.2f deg}', theta_deg), ...
     sprintf('\\rm{REQ DV    : %.4f km/s}', DeltaV), ...
     '\color[rgb]{1.0,0.4,0.0}\bf{ATTITUDE : NOMINAL}'}, ...
    'FontName', 'Consolas', 'FontSize', 9, 'Color', 'w', ...
    'BackgroundColor', [0.05 0.05 0.06 0.85], ... % Leggermente più trasparente
    'EdgeColor', [0.3 0.3 0.4], ...
    'Margin', 5, ...
    'Interpreter', 'tex', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle');

% --- IL TRUCCO MAGICO PER IL 3D ---
set(ax, 'SortMethod', 'depth');      % Forza MATLAB a calcolare chi sta davanti e chi dietro
set(t_hud, 'Layer', 'back');         % Impedisce al box di "saltare" sempre in primo piano
% --- 6e. MASTER HUD EVOLVED ---

% Colori HUD
c_cyan_hud  = [0.2, 1.0, 1.0];
c_white_hud = [0.9, 0.9, 0.9];


hud_str = {
    ' [ FLIGHT DYNAMICS SYSTEM v2.6 ]', ...
    ' ________________________________', ...
    ' ', ...
    sprintf('  > SEMIMAJOR AXIS :  %.2f km', a), ...
    sprintf('  > ECCENTRICITY   :  %.5f', e), ...
    sprintf('  > INITIAL INC    :  %.2f deg', rad2deg(i_i)), ...
    sprintf('  > FINAL INC      :  %.2f deg', rad2deg(i_f)), ...
    sprintf('  > ORIENTATION    :  %s', rot_str), ...
    ' ________________________________', ...
    ' ', ...
    sprintf('  > TARGET DV      :  %.4f km/s', DeltaV), ...
    sprintf('  > NODE ANGLE     :  %.2f deg', theta_deg), ...
    sprintf('  > INC DELTA      :  %.2f deg', rad2deg(abs(i_f-i_i))), ...
};

% Textbox Principale
h_hud = annotation('textbox', [0.02, 0.05, 0.22, 0.35], 'String', hud_str, ...
    'Color', c_white_hud, ...
    'EdgeColor', c_cyan_hud, ...
    'LineWidth', 1.5, ...
    'BackgroundColor', [0.02 0.02 0.03], ...
    'FaceAlpha', 0.85, ...
    'FontName', 'Consolas', ...
    'FontSize', 9, ...
    'Interpreter', 'none', ...
    'VerticalAlignment', 'middle');

% Effetto Cornice Doppia
annotation('rectangle', [0.018, 0.048, 0.224, 0.354], ...
    'EdgeColor', [c_cyan_hud 0.2], 'LineWidth', 0.5);

% Indicatore LIVE
annotation('textbox', [0.02, 0.405, 0.05, 0.02], 'String', 'o LIVE', ...
    'Color', [1 0.2 0.2], 'EdgeColor', 'none', ...
    'FontName', 'Consolas', 'FontSize', 8, 'FontWeight', 'bold');

% --- 6f. POST-PROCESSING ---
title('UPPER STAGE TRAJECTORY', 'Color', 'w', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
subtitle('HIGH FIDELITY ORBIT RENDER', 'Color', [0.5 0.5 0.5], 'FontName', 'Consolas', 'Interpreter', 'none');

% --- 6g. MASTER HUD LEGENDA ---
lgd = legend([h_orb_i, h_orb_f], ...
    'TextColor', 'w', ...
    'Color', [0.02 0.02 0.04 0.8], ... % Sfondo ultra-dark
    'EdgeColor', [0.3 0.3 0.5], ...
    'FontName', 'Consolas', ...
    'FontSize', 10, ...
    'Location', 'southoutside', ...   
    'Orientation', 'horizontal', ... % Dispone le voci una accanto all'altra
    'Box', 'on');

% Titolo della legenda 
lgd.Title.String = 'MISSION ASSETS IDENTIFIER';
lgd.Title.Color = [0.6 0.6 0.6];
lgd.Title.FontSize = 8;

grid off; 
drawnow;

% --- 7. AUTO-ZOOM, CAMERA SETUP & INTERACTIVITY UNLOCK ---
r_apocentro = a * (1 + e);
view_limit = r_apocentro * 1.15;

xlim([-view_limit, view_limit]);
ylim([-view_limit, view_limit]);
zlim([-view_limit, view_limit]);

set(ax, 'Projection', 'perspective');
axis vis3d; % impedisce all'orbita di "schiacciarsi" mentre ruoti

% SBLOCCO ROTAZIONE AUTOMATICA: Permette di usare il mouse anche durante il loop
h_rot = rotate3d(ax);
set(h_rot, 'Enable', 'on'); 

% Lighting
camlight('headlight');
lighting gouraud;
material dull; 
set(gcf, 'GraphicsSmoothing', 'on', 'InvertHardcopy', 'off');
