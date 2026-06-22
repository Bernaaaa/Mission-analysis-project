function scenery1_plot_standard(P1, P2, P3, P4, P5)
    % PLOT_STANDARD_SCENARIO - Grafico esplicito per la Strategia Standard
    
    mu = 398600.4418;
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
        surf(X_e*R_E, Y_e*R_E, Z_e*R_E, ...
             'CData', cdata, 'FaceColor', 'texturemap', 'EdgeColor', 'none', ...
             'FaceLighting', 'none', 'HandleVisibility', 'off');
        axis equal; 
    catch 
        surf(X_e*R_E, Y_e*R_E, Z_e*R_E, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'none', ...
             'FaceAlpha', 0.85, 'HandleVisibility', 'off');
        camlight('left'); lighting gouraud; material dull;
    end
     
    % --- PALETTE COLORI (Mappata sullo Scenario Standard) ---
    c_init  = [0.900 0.100 0.100]; % Rosso acceso (Init)
    c_t1    = [0.100 0.750 0.100]; % Verde smeraldo (Post-Plane)
    c_t2    = [0.000 0.750 0.900]; % Ciano brillante (Post-CP)
    c_bt    = [1.000 0.600 0.000]; % Arancio vivido (Bitangente)
    c_park  = [0.000 0.000 0.400]; % Blu navy scuro (Final Orbit)

    lw_arc = 0.8; 
    mk_sz  = 6;    
    
    vec_size = 1000;

    % Funzione inline per gestire i cavallotti sullo 0 (es. anomalia da 350° a 10°)
    get_th = @(t_i, t_f) linspace(t_i, t_f + (t_f < t_i)*2*pi, vec_size);

    % --- 2. PLOT DEGLI ARCHI REALI PERCORSI ---
    
    % Arco 1: Coasting su Orbita Iniziale
    th1 = get_th(P1.th_in, P1.th_in + 2*pi);
    R1 = zeros(3, vec_size); for j=1:vec_size, R1(:,j) = par2car(P1.a, P1.e, P1.i, P1.OM, P1.om, th1(j), mu); end
    plot3(R1(1,:), R1(2,:), R1(3,:), '-', 'Color', c_init, 'LineWidth', lw_arc, 'DisplayName', P1.name);
    
    % Arco 2: Trans Orbit (Post-Plane Change)
    th2 = get_th(P2.th_in, P2.th_out);
    R2 = zeros(3, vec_size); for j=1:vec_size, R2(:,j) = par2car(P2.a, P2.e, P2.i, P2.OM, P2.om, th2(j), mu); end
    plot3(R2(1,:), R2(2,:), R2(3,:), '-', 'Color', c_t1, 'LineWidth', lw_arc, 'DisplayName', P2.name);
    
    % Arco 3: Trans Orbit (Post-Pericenter Change)
    th3 = get_th(P3.th_in, P3.th_out);
    R3 = zeros(3, vec_size); for j=1:vec_size, R3(:,j) = par2car(P3.a, P3.e, P3.i, P3.OM, P3.om, th3(j), mu); end
    plot3(R3(1,:), R3(2,:), R3(3,:), '-', 'Color', c_t2, 'LineWidth', lw_arc, 'DisplayName', P3.name);
    
    % Arco 4: Trasferimento Bitangente
    th4 = get_th(P4.th_in, P4.th_out);
    R4 = zeros(3, vec_size); for j=1:vec_size, R4(:,j) = par2car(P4.a, P4.e, P4.i, P4.OM, P4.om, th4(j), mu); end
    plot3(R4(1,:), R4(2,:), R4(3,:), '-', 'Color', c_bt, 'LineWidth', lw_arc, 'DisplayName', P4.name);
    
    % Arco 5: Coasting su Orbita Finale
    th5 = get_th(P5.th_in, P5.th_out + 2*pi);
    R5 = zeros(3, vec_size); for j=1:vec_size, R5(:,j) = par2car(P5.a, P5.e, P5.i, P5.OM, P5.om, th5(j), mu); end
    plot3(R5(1,:), R5(2,:), R5(3,:), '-', 'Color', c_park, 'LineWidth', lw_arc, 'DisplayName', P5.name);

    % --- 3. PUNTI DI MANOVRA ---
    mkOpt = {'o', 'LineWidth', 1.5, 'MarkerSize', mk_sz};
    
    plot3(R1(1,1), R1(2,1), R1(3,1), mkOpt{:}, 'MarkerEdgeColor', c_init, 'DisplayName', P1.maneuver);
    plot3(R2(1,1), R2(2,1), R2(3,1), mkOpt{:}, 'MarkerEdgeColor', c_t1, 'DisplayName', P2.maneuver);
    plot3(R3(1,1), R3(2,1), R3(3,1), mkOpt{:}, 'MarkerEdgeColor', c_t2, 'DisplayName', P3.maneuver);
    plot3(R4(1,1), R4(2,1), R4(3,1), mkOpt{:}, 'MarkerEdgeColor', c_bt, 'DisplayName', P4.maneuver);
    plot3(R5(1,1), R5(2,1), R5(3,1), mkOpt{:}, 'MarkerEdgeColor', c_park, 'DisplayName', P5.maneuver);
    
    % Target Finale
    plot3(R5(1,end), R5(2,end), R5(3,end), 'o', 'MarkerEdgeColor', 'k', ...
          'LineWidth', 1.5, 'MarkerSize', mk_sz, 'DisplayName', 'Final Target');

    % --- 4. DETTAGLI EDITORIALI E ZOOM OTTIMIZZATO ---
    axis equal;
    
    % Calcola lo zoom unicamente sull'orbita più ampia
    max_r = max([P1.a*(1+P1.e), P2.a*(1+P2.e), P3.a*(1+P3.e), P4.a*(1+P4.e), P5.a*(1+P5.e)]);
    lim = max_r * 1.1; 
    
    xlim([-lim, lim]); ylim([-lim, lim]); zlim([-lim, lim]);
    
    xlabel('$X$ [km]', 'FontSize', 12, 'Interpreter', 'latex');
    ylabel('$Y$ [km]', 'FontSize', 12, 'Interpreter', 'latex');
    zlabel('$Z$ [km]', 'FontSize', 12, 'Interpreter', 'latex');
    title('\textbf{Scenario 1: Standard Transfer}', 'FontSize', 14, 'Interpreter', 'latex');
    
    view(145, 20);
    
    lgd = legend('show', 'Location', 'best');
    lgd.FontSize = 10;
    lgd.Title.String = 'Mission Phases \& Maneuvers';
    lgd.Title.Interpreter = 'latex';
    lgd.Box = 'on';

end