function scenery1_plot(GTO, Park, Trans1, Trans2)
    
    mu = GTO.mu;
    R_E = 6378.137;   % [km] Raggio terrestre
    
    % Configurazione globale LaTeX
    set(groot, 'defaulttextinterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    
    figure('Color', 'w', 'WindowState', 'maximized');
    hold on; grid on; box on;
    ax = gca; 
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.05;
    
    % --- 1. TERRA 3D CON TEXTURE ---
    [X_e, Y_e, Z_e] = sphere(100);
    
    try 
        cdata = imread('earth_texture.png'); 
        cdata = flipud(cdata); 
        
        % Disegno della Terra
        surf(X_e*R_E, Y_e*R_E, Z_e*R_E, ...
             'CData', cdata, ...
             'FaceColor', 'texturemap', ...
             'EdgeColor', 'none', ...
             'FaceLighting', 'none', ...
             'HandleVisibility', 'off');
        
        axis equal; 
    catch 
        warning('Texture per la Terra non trovata. Disegno una sfera semplice.');
        surf(X_e*R_E, Y_e*R_E, Z_e*R_E, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'none', ...
             'FaceAlpha', 0.85, 'HandleVisibility', 'off');
        camlight('left');
        lighting gouraud;
        material dull;
    end
     
    % --- PALETTE COLORI ---
    c_gto   = [0.900 0.100 0.100]; % Rosso acceso
    c_bt1   = [1.000 0.600 0.000]; % Arancio vivido
    c_t1    = [0.100 0.750 0.100]; % Verde smeraldo
    c_mid   = [0.000 0.750 0.900]; % Ciano brillante
    c_t2    = [0.750 0.100 0.850]; % Magenta/Viola forte
    c_bt2   = [0.100 0.300 1.000]; % Blu elettrico
    c_park  = [0.000 0.000 0.400]; % Blu navy scuro 

    % Spessori
    lw_arc = 0.8;  % Linee sottilissime
    mk_sz  = 6;    % Dimensione marker invariata
    
    % Recupero om_mid
    if isfield(Trans1, 'om_mid'); om_mid = Trans1.om_mid; else; om_mid = Trans1.om; end
    
    % Calcolo rami bitangenti
    rp_GTO = GTO.a * (1 - GTO.e); ra_T1 = Trans1.a * (1 + Trans1.e);
    a_bt1 = (rp_GTO + ra_T1) / 2; e_bt1 = (ra_T1 - rp_GTO) / (ra_T1 + rp_GTO);
    
    ra_T2 = Trans2.a * (1 + Trans2.e); rp_Park = Park.a * (1 - Park.e);
    a_bt2 = (rp_Park + ra_T2) / 2; e_bt2 = (ra_T2 - rp_Park) / (ra_T2 + rp_Park);

    % --- 2. PLOT DEGLI ARCHI REALI PERCORSI ---
    
    % Arco 1: Coasting su GTO
    th_gto = linspace(0, 2*pi, 100);
    R_gto_arc = zeros(3, 100); for j=1:100, R_gto_arc(:,j) = par2car(GTO.a, GTO.e, GTO.i, GTO.OM, GTO.om, th_gto(j), mu); end
    plot3(R_gto_arc(1,:), R_gto_arc(2,:), R_gto_arc(3,:), '-', 'Color', c_gto, 'LineWidth', lw_arc, 'DisplayName', 'GTO Coasting');
    
    % Arco 2: Primo Ramo Bitangente
    th_bt1 = linspace(0, pi, 100);
    R_bt1_arc = zeros(3, 100); for j=1:100, R_bt1_arc(:,j) = par2car(a_bt1, e_bt1, GTO.i, GTO.OM, GTO.om, th_bt1(j), mu); end
    plot3(R_bt1_arc(1,:), R_bt1_arc(2,:), R_bt1_arc(3,:), '-', 'Color', c_bt1, 'LineWidth', lw_arc, 'DisplayName', 'BT Leg 1');
    
    % Estrazione anomalie manovre
    th_plane = pi; if isfield(Trans1, 'theta_plane'), th_plane = Trans1.theta_plane; end
    th_omega_i = pi; if isfield(Trans1, 'thi_omega'), th_omega_i = Trans1.thi_omega; end
    th_omega_f = pi; if isfield(Trans2, 'thf_omega'), th_omega_f = Trans2.thf_omega; end
    
    % Arco 3: Su Trans1
    th_t1_arc = linspace(pi, th_plane, 80);
    R_t1_arc = zeros(3, length(th_t1_arc)); for j=1:length(th_t1_arc), R_t1_arc(:,j) = par2car(Trans1.a, Trans1.e, Trans1.i, Trans1.OM, Trans1.om, th_t1_arc(j), mu); end
    plot3(R_t1_arc(1,:), R_t1_arc(2,:), R_t1_arc(3,:), '-', 'Color', c_t1, 'LineWidth', lw_arc, 'DisplayName', 'Trans Orbit (Pre-Plane)');
    
    % Arco 4: Tratto intermedio
    th_mid_arc = linspace(th_plane, th_omega_i, 80);
    R_mid_arc = zeros(3, length(th_mid_arc)); for j=1:length(th_mid_arc), R_mid_arc(:,j) = par2car(Trans1.a, Trans1.e, Park.i, Park.OM, om_mid, th_mid_arc(j), mu); end
    plot3(R_mid_arc(1,:), R_mid_arc(2,:), R_mid_arc(3,:), '-', 'Color', c_mid, 'LineWidth', lw_arc, 'DisplayName', 'Trans Orbit (Pre-$\omega$)');
    
    % Arco 5: Su Trans2
    th_t2_arc = linspace(th_omega_f, pi, 80);
    R_t2_arc = zeros(3, length(th_t2_arc)); for j=1:length(th_t2_arc), R_t2_arc(:,j) = par2car(Trans2.a, Trans2.e, Trans2.i, Trans2.OM, Trans2.om, th_t2_arc(j), mu); end
    plot3(R_t2_arc(1,:), R_t2_arc(2,:), R_t2_arc(3,:), '-', 'Color', c_t2, 'LineWidth', lw_arc, 'DisplayName', 'Trans Orbit (Final)');
    
    % Arco 6: Secondo Ramo Bitangente
    th_bt2 = linspace(pi, 2*pi, 100);
    R_bt2_arc = zeros(3, 100); for j=1:100, R_bt2_arc(:,j) = par2car(a_bt2, e_bt2, Park.i, Park.OM, Park.om, th_bt2(j), mu); end
    plot3(R_bt2_arc(1,:), R_bt2_arc(2,:), R_bt2_arc(3,:), '-', 'Color', c_bt2, 'LineWidth', lw_arc, 'DisplayName', 'BT Leg 2');
    
    % Arco 7: Coasting su Parking Orbit
    th_park = linspace(0, Park.th, 100);
    R_park_arc = zeros(3, 100); for j=1:100, R_park_arc(:,j) = par2car(Park.a, Park.e, Park.i, Park.OM, Park.om, th_park(j), mu); end
    plot3(R_park_arc(1,:), R_park_arc(2,:), R_park_arc(3,:), '-', 'Color', c_park, 'LineWidth', lw_arc, 'DisplayName', 'Park Coasting');

    % --- 3. PUNTI DI MANOVRA ---
    mkOpt = {'o', 'LineWidth', 1.5, 'MarkerSize', mk_sz};
    
    r_start = par2car(GTO.a, GTO.e, GTO.i, GTO.OM, GTO.om, GTO.th, mu);
    plot3(r_start(1), r_start(2), r_start(3), mkOpt{:}, 'MarkerEdgeColor', c_gto, 'DisplayName', 'Initial Pos');
    
    r_m1 = par2car(GTO.a, GTO.e, GTO.i, GTO.OM, GTO.om, 0, mu);
    plot3(r_m1(1), r_m1(2), r_m1(3), mkOpt{:}, 'MarkerEdgeColor', c_gto, 'DisplayName', '$\Delta V_{1BT}$');
    
    r_m2 = par2car(a_bt1, e_bt1, GTO.i, GTO.OM, GTO.om, pi, mu);
    plot3(r_m2(1), r_m2(2), r_m2(3), mkOpt{:}, 'MarkerEdgeColor', c_bt1, 'DisplayName', '$\Delta V_{2BT}$');
    
    r_mP = par2car(Trans1.a, Trans1.e, Trans1.i, Trans1.OM, Trans1.om, th_plane, mu);
    plot3(r_mP(1), r_mP(2), r_mP(3), mkOpt{:}, 'MarkerEdgeColor', c_t1, 'DisplayName', '$\Delta V_{plane}$');
    
    r_mOm = par2car(Trans1.a, Trans1.e, Park.i, Park.OM, om_mid, th_omega_i, mu);
    plot3(r_mOm(1), r_mOm(2), r_mOm(3), mkOpt{:}, 'MarkerEdgeColor', c_mid, 'DisplayName', '$\Delta V_{\omega}$');
    
    r_m3 = par2car(Trans2.a, Trans2.e, Trans2.i, Trans2.OM, Trans2.om, pi, mu);
    plot3(r_m3(1), r_m3(2), r_m3(3), mkOpt{:}, 'MarkerEdgeColor', c_t2, 'DisplayName', '$\Delta V_{3BT}$');
    
    r_m4 = par2car(a_bt2, e_bt2, Park.i, Park.OM, Park.om, 2*pi, mu);
    plot3(r_m4(1), r_m4(2), r_m4(3), mkOpt{:}, 'MarkerEdgeColor', c_bt2, 'DisplayName', '$\Delta V_{4BT}$');
    
    r_end = par2car(Park.a, Park.e, Park.i, Park.OM, Park.om, Park.th, mu);
    plot3(r_end(1), r_end(2), r_end(3), mkOpt{:}, 'MarkerEdgeColor', c_park, 'DisplayName', 'Final Pos');

    % --- 4. DETTAGLI EDITORIALI ---
    max_ap = max([GTO.a*(1+GTO.e), Trans1.a*(1+Trans1.e), Park.a*(1+Park.e), a_bt1*(1+e_bt1), a_bt2*(1+e_bt2)]);
    axis equal;
    lim = max_ap * 1.05; 
    xlim([-lim, lim]); ylim([-lim, lim]); zlim([-lim, lim]);
    
    xlabel('$X$ [km]', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('$Y$ [km]', 'FontSize', 12, 'Interpreter', 'latex');
    zlabel('$Z$ [km]', 'FontSize', 12, 'Interpreter', 'latex');
    title('\textbf{Scenario 1: Bielliptic Transfer \& Nodal Maneuvers}', 'FontSize', 14, 'Interpreter', 'latex');
    
    view(145, 20);
    
    lgd = legend('show', 'Location', 'best');
    lgd.FontSize = 10;
    lgd.Title.String = 'Mission Phases \& Maneuvers';
    lgd.Title.Interpreter = 'latex';
    lgd.Box = 'on';

    % --- 5. RIQUADRO ZOOM (INSET PLOT) SULLA TERRA ---
    
    % Posizione esatta in basso a sinistra
    ax_zoom = axes('Position', [-0.05 0.25 0.4 0.4]); 
    
    % Clona gli elementi dal grafico principale
    copyobj(allchild(ax), ax_zoom);
    
    % Copia l'angolo di visuale
    view(ax_zoom, 145, 20); 
    grid(ax_zoom, 'on'); 
    box(ax_zoom, 'on');
    
    % --- IL FIX PER LA TERRA SCHIACCIATA ---
    axis(ax_zoom, 'equal'); 
    
    % Calcola i limiti di zoom 
    lim_zoom = max(GTO.a*(1+GTO.e), Park.a*(1+Park.e)) * 1.15;
    
    % Applica i nuovi limiti
    xlim(ax_zoom, [-lim_zoom, lim_zoom]);
    ylim(ax_zoom, [-lim_zoom, lim_zoom]);
    zlim(ax_zoom, [-lim_zoom, lim_zoom]);
    
    % --- IL FIX PER LA SOVRAPPOSIZIONE ---
    set(ax_zoom, 'Color', 'w'); % Sfondo bianco solido
    
    % Estetica del riquadro
    title(ax_zoom, '\textbf{Zoom: Low Earth Proximity}', 'Interpreter', 'latex', 'FontSize', 12);
    xticklabels(ax_zoom, {}); yticklabels(ax_zoom, {}); zticklabels(ax_zoom, {});
    ax_zoom.GridAlpha = 0.2;

end