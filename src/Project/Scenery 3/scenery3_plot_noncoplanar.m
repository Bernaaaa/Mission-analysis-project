function scenery3_plot_noncoplanar(Park, Hyp_peri, Hyp_opt, mu, R_earth)
    
    % Configurazione globale LaTeX
    set(groot, 'defaulttextinterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    % Palette colori
    c_park  = [0.000 0.000 0.400]; % Blu (Orbita Parcheggio)
    c_hyp   = [0.900 0.100 0.100]; % Rosso (Iperbole)
    c_man   = [1.000 0.600 0.000]; % Arancio (Punto Manovra)
    
    lw_arc = 1.5; 
    lw_hyp = 2.0;
    mk_sz  = 7;
    vec_size = 1000; % Risoluzione archi

    function draw_scenario(Hyp_data, title_str)
        figure('Color', 'w', 'WindowState', 'maximized', 'Name', title_str);
        hold on; grid on; box on; axis equal;
        ax = gca; 
        ax.GridAlpha = 0.15; 
        ax.MinorGridAlpha = 0.05;

        % 1. Terra 3D
        [X_e, Y_e, Z_e] = sphere(100);
        try 
            cdata = imread('earth_texture.png'); 
            cdata = flipud(cdata); 
            surf(X_e*R_earth, Y_e*R_earth, Z_e*R_earth, 'CData', cdata, 'FaceColor', 'texturemap', ...
                 'EdgeColor', 'none', 'FaceLighting', 'none', 'HandleVisibility', 'off');
        catch 
            surf(X_e*R_earth, Y_e*R_earth, Z_e*R_earth, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'none', ...
                 'FaceAlpha', 0.85, 'HandleVisibility', 'off');
            camlight('left'); lighting gouraud; material dull;
        end

        % 2. Orbita di Parcheggio
        th_park = linspace(0, 2*pi, vec_size);
        R_park = zeros(3, vec_size);
        for j=1:vec_size
            R_park(:,j) = par2car(Park.a, Park.e, Park.i, Park.OM, Park.om, th_park(j), mu); 
        end
        plot3(R_park(1,:), R_park(2,:), R_park(3,:), '-', 'Color', c_park, 'LineWidth', lw_arc, 'DisplayName', 'Park Orbit');

        % 3. Iperbole di Fuga
        th_inf_H = acos(-1 / Hyp_data.e);
        th_hyp = linspace(-th_inf_H*0.85, th_inf_H*0.85, vec_size); 
        R_hyp = zeros(3, vec_size);
        for j=1:vec_size
            R_hyp(:,j) = par2car(Hyp_data.a, Hyp_data.e, Hyp_data.i, Hyp_data.OM, Hyp_data.om, th_hyp(j), mu); 
        end
        plot3(R_hyp(1,:), R_hyp(2,:), R_hyp(3,:), '-', 'Color', c_hyp, 'LineWidth', lw_hyp, 'DisplayName', 'Hyperbola of Escape');

        % 4. Punto di Iniezione (Manovra)
        R_inj = par2car(Park.a, Park.e, Park.i, Park.OM, Park.om, Hyp_data.th_in, mu);
        plot3(R_inj(1), R_inj(2), R_inj(3), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c_man, ...
              'LineWidth', 1.5, 'MarkerSize', mk_sz, 'DisplayName', 'Injection Point ($\Delta v$)');

              
        lim_E = max(vecnorm(R_park)) * 1.5; 
        xlim([-lim_E, lim_E]); ylim([-lim_E, lim_E]); zlim([-lim_E, lim_E]);
        xlabel('$X$ [km]', 'FontSize', 12); ylabel('$Y$ [km]', 'FontSize', 12); zlabel('$Z$ [km]', 'FontSize', 12);
        title(['\textbf{', title_str, '}'], 'FontSize', 18);
        
        view(110, 20); 
        
        lgd = legend('show', 'Location', 'northeast'); 
        lgd.Box = 'on';
    end

    draw_scenario(Hyp_peri, 'Phase 1: Earth Escape (Pericenter Injection)');
    draw_scenario(Hyp_opt,  'Phase 1: Earth Escape (Optimal Injection)');
end