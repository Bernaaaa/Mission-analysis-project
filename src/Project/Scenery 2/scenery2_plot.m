function scenery2_plot(Earth, Nyx, Trans, Const)
    % --- CONFIGURAZIONE COLORI ---
    palette.terra = [0, 0.45, 0.74]; 
    palette.nyx   = [0.85, 0.33, 0.1]; 
    palette.trans = [0.7, 0, 1]; 

    thi_t = Const.th_i_t;
    if thi_t > Const.th_f_t
        thi_t = thi_t - 2*pi;
    end
    
    % Generazione coordinate (AU)
    R_E = get_coords(Earth, linspace(0, 2*pi, 800), Const.AU);
    R_N = get_coords(Nyx,   linspace(0, 2*pi, 800), Const.AU);
    R_T_full = get_coords(Trans, linspace(0, 2*pi, 800), Const.AU);
    R_T_arc  = get_coords(Trans, linspace(thi_t, Const.th_f_t, 400), Const.AU);

    % --- PLOT 1: GEOMETRIA COMPLETA ---
    start_fig('Mission Geometry');
    draw_sun_textured();
    plot3(R_E(1,:), R_E(2,:), R_E(3,:), 'Color', palette.terra, 'LineWidth', 1, 'DisplayName', 'Earth');
    plot3(R_N(1,:), R_N(2,:), R_N(3,:), 'Color', palette.nyx,   'LineWidth', 1, 'DisplayName', 'Nyx');
    plot3(R_T_arc(1,:), R_T_arc(2,:), R_T_arc(3,:), 'Color', palette.trans, 'LineWidth', 3, 'DisplayName', 'Transfer Arc');
    add_markers(Const, palette);
    apply_style(gca, 3.2, 'Earth to 3908 Nyx Mission Geometry');

    % --- PLOT 2: TERRA ---
    start_fig('Solo Earth');
    draw_sun_textured();
    plot3(R_E(1,:), R_E(2,:), R_E(3,:), 'Color', palette.terra, 'LineWidth', 1.5);
    apply_style(gca, 1.5, 'Earth Heliocentric Orbit');

    % --- PLOT 3: NYX ---
    start_fig('Solo Nyx');
    draw_sun_textured();
    plot3(R_N(1,:), R_N(2,:), R_N(3,:), 'Color', palette.nyx, 'LineWidth', 1.5);
    apply_style(gca, 3.2, '3908 Nyx Heliocentric Orbit');

    % --- PLOT 4: TRASFERIMENTO ---
    start_fig('Solo Transfer');
    draw_sun_textured();
    plot3(R_T_full(1,:), R_T_full(2,:), R_T_full(3,:), 'Color', palette.trans, 'LineWidth', 1.5);
    plot3(R_T_arc(1,:), R_T_arc(2,:), R_T_arc(3,:), 'Color', palette.trans, 'LineWidth', 1.5);
    add_markers(Const, palette);
    apply_style(gca, 2.5, 'Transfer Orbit Geometry');
end

% --- FUNZIONI DI SUPPORTO ---

function draw_sun_textured()
    [x, y, z] = sphere(25); 
    r_sun = 0.13; % Dimensione visibile in AU
    
    try
        % Lettura texture
        img = imread('sun_texture.jpg');
        s = surface(x*r_sun, y*r_sun, z*r_sun);
        
        % Applicazione texture "Lite"
        set(s, 'FaceColor', 'texturemap', ...
               'CData', img, ...
               'EdgeColor', 'none', ...
               'FaceLighting', 'none', ...
               'HandleVisibility', 'off');
    catch
        % se manca il file: sfera gialla semplice
        surface(x*r_sun, y*r_sun, z*r_sun, 'FaceColor', [1 0.8 0], 'EdgeColor', 'none');
    end
end

function R = get_coords(Orbit, th_vec, AU)
    R = zeros(3, length(th_vec));
    for j = 1:length(th_vec)
        [r, ~] = par2car(Orbit.a, Orbit.e, Orbit.i, Orbit.OM, Orbit.om, th_vec(j), Orbit.mu);
        R(:,j) = r / AU;
    end
end

function f = start_fig(name)
    f = figure('Color', 'w', 'Name', name);
    set(f, 'Units', 'normalized', 'Position', [0.2 0.2 0.4 0.5]);
    hold on; axis equal; grid on; view(135, 30);
end

function add_markers(C, p)
    plot3(C.r1_opt(1)/C.AU, C.r1_opt(2)/C.AU, C.r1_opt(3)/C.AU, 'o', 'MarkerFaceColor', p.terra, 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off');
    plot3(C.r2_opt(1)/C.AU, C.r2_opt(2)/C.AU, C.r2_opt(3)/C.AU, 'd', 'MarkerFaceColor', p.nyx, 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off');
end

function apply_style(ax, lim, tit)
    ax.TickLabelInterpreter = 'latex';
    ax.GridAlpha = 0.1;
    xlabel('$X$ [AU]', 'Interpreter', 'latex');
    ylabel('$Y$ [AU]', 'Interpreter', 'latex');
    zlabel('$Z$ [AU]', 'Interpreter', 'latex');
    title(['\textbf{', tit, '}'], 'Interpreter', 'latex', 'FontSize', 12);
    xlim([-lim lim]); ylim([-lim lim]); zlim([-1 1]);
end