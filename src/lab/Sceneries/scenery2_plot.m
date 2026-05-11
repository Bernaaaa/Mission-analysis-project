function scenery2_plot(Earth, Nyx, Trans, Const)
    % Scenery2_plot: Genera i 4 plot (Global, Terra, Nyx, Trasferimento)
    % 
    % INPUT: 
    %   Earth, Nyx, Trans: struct con i parametri [a, e, i, OM, om, mu]
    %   Const: struct con [AU, r1_opt, r2_opt, th_i_t, th_f_t]

    % --- CONFIGURAZIONE COLORI E STILE ---
    palette.terra = [0, 0.45, 0.74]; 
    palette.nyx   = [0.85, 0.33, 0.1]; 
    palette.trans = [0.47, 0.67, 0.19]; 
    palette.sun   = [0.93, 0.69, 0.13];

    % Generazione coordinate (AU)
    R_E = get_coords(Earth, linspace(0, 2*pi, 1000), Const.AU);
    R_N = get_coords(Nyx,   linspace(0, 2*pi, 1000), Const.AU);
    R_T_full = get_coords(Trans, linspace(0, 2*pi, 1000), Const.AU);
    R_T_arc  = get_coords(Trans, linspace(Const.th_i_t, Const.th_f_t, 500), Const.AU);

    % --- PLOT 1: GEOMETRIA COMPLETA ---
    draw_sun(palette.sun);
    plot3(R_E(1,:), R_E(2,:), R_E(3,:), 'Color', palette.terra, 'LineWidth', 1, 'DisplayName', 'Earth');
    plot3(R_N(1,:), R_N(2,:), R_N(3,:), 'Color', palette.nyx,   'LineWidth', 1, 'DisplayName', 'Nyx');
    plot3(R_T_arc(1,:), R_T_arc(2,:), R_T_arc(3,:), 'Color', palette.trans, 'LineWidth', 3, 'DisplayName', 'Transfer Arc');
    add_markers(Const, palette);
    apply_style(gca, 3.2, 'Earth to 3908 Nyx Mission Geometry');

    % --- PLOT 2: TERRA ---
    start_fig('Solo Earth');
    draw_sun(palette.sun);
    plot3(R_E(1,:), R_E(2,:), R_E(3,:), 'Color', palette.terra, 'LineWidth', 1.5);
    apply_style(gca, 1.5, 'Earth Heliocentric Orbit');

    % --- PLOT 3: NYX ---
    start_fig('Solo Nyx');
    draw_sun(palette.sun);
    plot3(R_N(1,:), R_N(2,:), R_N(3,:), 'Color', palette.nyx, 'LineWidth', 1.5);
    apply_style(gca, 3.2, '3908 Nyx Heliocentric Orbit');

    % --- PLOT 4: TRASFERIMENTO ---
    start_fig('Solo Transfer');
    draw_sun(palette.sun);
    plot3(R_T_full(1,:), R_T_full(2,:), R_T_full(3,:), 'Color', palette.trans, 'LineWidth', 1.5, 'DisplayName', 'Full Ellipse');
    add_markers(Const, palette);
    apply_style(gca, 2.5, 'Transfer Orbit Geometry');
end

% --- FUNZIONI DI SUPPORTO (LOCALI) ---

function R = get_coords(Orbit, th_vec, AU)
    R = zeros(3, length(th_vec));
    for j = 1:length(th_vec)
        [r, ~] = par2car(Orbit.a, Orbit.e, Orbit.i, Orbit.OM, Orbit.om, th_vec(j), Orbit.mu);
        R(:,j) = r / AU;
    end
end

function f = start_fig(name)
    f = figure('Color', 'w', 'Name', name);
    set(f, 'Units', 'normalized', 'Position', [0.2 0.2 0.5 0.6]);
    hold on; axis equal; grid on; view(135, 30);
end

function draw_sun(c)
    plot3(0,0,0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

function add_markers(C, p)
    % Partenza P1
    plot3(C.r1_opt(1)/C.AU, C.r1_opt(2)/C.AU, C.r1_opt(3)/C.AU, 'ok', 'MarkerFaceColor', p.terra, 'MarkerSize', 8, 'HandleVisibility', 'off');
    % Arrivo P2
    plot3(C.r2_opt(1)/C.AU, C.r2_opt(2)/C.AU, C.r2_opt(3)/C.AU, 'sk', 'MarkerFaceColor', p.nyx, 'MarkerSize', 8, 'HandleVisibility', 'off');
end

function apply_style(ax, lim, tit)
    ax.TickLabelInterpreter = 'latex';
    ax.GridAlpha = 0.15;
    xlabel('$X$ [AU]', 'Interpreter', 'latex');
    ylabel('$Y$ [AU]', 'Interpreter', 'latex');
    zlabel('$Z$ [AU]', 'Interpreter', 'latex');
    title(['\textbf{', tit, '}'], 'Interpreter', 'latex', 'FontSize', 14);
    xlim([-lim lim]); ylim([-lim lim]); zlim([-lim/3 lim/3]);
end