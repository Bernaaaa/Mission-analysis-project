function scenery3_plot_coplanar(E, N)
    
    % Configurazione globale LaTeX
    set(groot, 'defaulttextinterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    % Palette colori
    c_park  = [0.000 0.000 0.400]; % Blu (Orbita Parcheggio)
    c_hyp   = [0.900 0.100 0.100]; % Rosso (Iperboli)
    c_circ  = [0.100 0.750 0.100]; % Verde (Cattura circolare)
    c_man   = [1.000 0.600 0.000]; % Arancio (Punti Manovra)
    
    lw_arc = 1.5; 
    lw_hyp = 2.0;
    mk_sz  = 7;
    vec_size = 750; % Risoluzione archi

    % FIGURA 1: EARTH DEPARTURE
    figure('Color', 'w', 'WindowState', 'maximized', 'Name', 'Phase 1: Earth Escape');
    hold on; grid on; box on; axis equal;
    ax1 = gca; 
    ax1.GridAlpha = 0.15; 
    ax1.MinorGridAlpha = 0.05;

    % 1. Terra 3D
    [X_e, Y_e, Z_e] = sphere(100);
    try 
        cdata = imread('earth_texture.png'); 
        cdata = flipud(cdata); 
        surf(X_e*E.R, Y_e*E.R, Z_e*E.R, 'CData', cdata, 'FaceColor', 'texturemap', ...
             'EdgeColor', 'none', 'FaceLighting', 'none', 'HandleVisibility', 'off');
    catch 
        surf(X_e*E.R, Y_e*E.R, Z_e*E.R, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'none', ...
             'FaceAlpha', 0.85, 'HandleVisibility', 'off');
        camlight('left'); lighting gouraud; material dull;
    end

    % 2. Orbita di Parcheggio
    th_park = linspace(0, 2*pi, vec_size);
    R_park = zeros(3, vec_size);

    for j=1:vec_size
        R_park(:,j) = par2car(E.a_park, E.e_park, E.i, E.OM, E.om, th_park(j), E.mu); 
    end

    plot3(R_park(1,:), R_park(2,:), R_park(3,:), '-', 'Color', c_park, 'LineWidth', lw_arc, 'DisplayName', 'Parking Orbit');

    % 3. Iperbole di Fuga
    th_inf_E = acos(-1 / E.e_hyp);
    th_hyp = linspace(-th_inf_E*0.85, th_inf_E*0.85, vec_size); 
    R_hyp = zeros(3, vec_size);
    R_departure = par2car(E.a_hyp, E.e_hyp, E.i, E.OM, E.om, 0, E.mu);

    for j=1:vec_size
        R_hyp(:,j) = par2car(E.a_hyp, E.e_hyp, E.i, E.OM, E.om, th_hyp(j), E.mu); 
    end

    plot3(R_hyp(1,:), R_hyp(2,:), R_hyp(3,:), '-', 'Color', c_hyp, 'LineWidth', lw_hyp, ...
     'DisplayName', 'Departure Hyperbola');

    % 4. Manovra (Pericentro)
    plot3(R_departure(1), R_departure(2), R_departure(3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_man, ...
          'LineWidth', 1.5, 'MarkerSize', mk_sz, 'DisplayName', '$\Delta V$ Departure');

    lim_E = max(vecnorm(R_park)) * 1.1; 
    xlim([-lim_E, lim_E]); ylim([-lim_E, lim_E]); zlim([-lim_E, lim_E]);
    xlabel('$X$ [km]', 'FontSize', 12); ylabel('$Y$ [km]', 'FontSize', 12); zlabel('$Z$ [km]', 'FontSize', 12);
    title('\textbf{Phase 1: Earth Escape}', 'FontSize', 18);
    view(145, 0); % Vista laterale
    lgd1 = legend('show', 'Location', 'northeast'); 
    lgd1.Box = 'on';

    % FIGURA 2: NYX ARRIVAL
    figure('Color', 'w', 'WindowState', 'maximized', 'Name', 'Phase 3: Nyx Capture');
    hold on; grid on; box on; axis equal;
    ax2 = gca; 
    ax2.GridAlpha = 0.15; 
    ax2.MinorGridAlpha = 0.05;

    % Asteroide
    radius = N.R;
    [X_n, Y_n, Z_n] = sphere(100);

    try 
        cdata = imread('asteroid_texture.jpg'); 
        cdata = flipud(cdata); 
        surf(X_n*radius, Y_n*radius, Z_n*radius, 'CData', cdata, 'FaceColor', 'texturemap', ...
             'EdgeColor', 'none', 'FaceLighting', 'none', 'HandleVisibility', 'off');
    catch 
        surf(X_n*radius, Y_n*radius, Z_n*radius, 'FaceColor', [0.5 0.5 0.5], ...
         'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'DisplayName', 'Nyx (Scaled)');
        camlight('left'); material dull;
    end

    % Parametri orbitali locali
    i_n = 0; 
    OM_n = 0; 
    om_n = 0;

    % 2. Orbita Circolare di Cattura
    th_circ = linspace(0, 2*pi, vec_size);
    R_circ = zeros(3, vec_size);

    for j=1:vec_size
        R_circ(:,j) = par2car(N.r_circ, 0, i_n, OM_n, om_n, th_circ(j), N.mu); 
    end

    plot3(R_circ(1,:), R_circ(2,:), R_circ(3,:), '-', 'Color', c_circ, 'LineWidth', lw_arc, 'DisplayName', 'Capture Orbit');

    % 3. Iperbole di Arrivo 
    th_inf_N = acos(-1 / N.e_hyp);
    th_hyp_n = linspace(-th_inf_N*0.85, th_inf_N*0.85, vec_size);
    R_hyp_n = zeros(3, vec_size);
    R_arrival = par2car(N.a_hyp, N.e_hyp, i_n, OM_n, om_n, 0, N.mu);

    for j=1:vec_size 
        R_hyp_n(:,j) = par2car(N.a_hyp, N.e_hyp, i_n, OM_n, om_n, th_hyp_n(j), N.mu); 
    end
    
    plot3(R_hyp_n(1,:), R_hyp_n(2,:), R_hyp_n(3,:), '-', 'Color', c_hyp, 'LineWidth', lw_hyp, ...
     'DisplayName', 'Arrival Hyperbola');

    % 4. Manovra (Pericentro di arrivo)
    plot3(R_arrival(1), R_arrival(2), R_arrival(3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_man, ...
          'LineWidth', 1.5, 'MarkerSize', mk_sz, 'DisplayName', '$\Delta V$ Capture');

    lim_N = max(vecnorm(R_circ)) * 2.0; 
    xlim([-lim_N, lim_N]); ylim([-lim_N, lim_N]); zlim([-lim_N, lim_N]);
    xlabel('$X$ [km]', 'FontSize', 12); ylabel('$Y$ [km]', 'FontSize', 12); zlabel('$Z$ [km]', 'FontSize', 12);
    title('\textbf{Phase 3: Nyx Capture}', 'FontSize', 18);
    view(45, 0);
    lgd2 = legend('show', 'Location', 'northeast'); 
    lgd2.Box = 'on';

end