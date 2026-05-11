function [best_th_earth, best_th_nyx, best_om_tn, best_dv] = grid_search_s2(a_i, e_i, i_i, OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu, R_sun, grid_size)

    th_i = linspace(0, 2*pi, grid_size); % true anomaly of the initial orbit varying from 0 to 360 degrees
    th_f = linspace(0, 2*pi, grid_size); % true anomaly of the final orbit varying from 0 to 360 degrees
    om_tn = linspace(0, 2*pi, grid_size); % pericenter argument of the transfer orbit varying from 0 to 360 degrees

    r1 = zeros(3, grid_size);
    v_i = zeros(3, grid_size);
    r2 = zeros(3, grid_size);
    v_f = zeros(3, grid_size);

    r1_n = zeros(1, grid_size);
    r2_n = zeros(1, grid_size);

    best_dv = inf;

    for i = 1:grid_size
        [r1(:,i), v_i(:,i)] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i(i), mu);
        [r2(:,i), v_f(:,i)] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f(i), mu); 
        r1_n(i) = norm(r1(:,i));
        r2_n(i) = norm(r2(:,i));
    end

    for j = 1:grid_size

        r1_current = r1(:,j);  v_i_current = v_i(:,j);  r1_n_current = r1_n(j);

        for k = 1:grid_size

            r2_current = r2(:,k);  v_f_current = v_f(:,k);  r2_n_current = r2_n(k);
                            
            p_vect_r1r2 = cross(r1_current, r2_current);
            h_t = p_vect_r1r2 / norm(p_vect_r1r2);
            i_t = acos(h_t(3));

            p_vect_kht = cross([0; 0; 1], p_vect_r1r2);
            N_t = p_vect_kht / norm(p_vect_kht);

            if N_t(2) >= 0
                OM_t = acos(N_t(1));
            else
                OM_t = 2*pi - acos(N_t(1));
            end

            R3_OM = [ cos(OM_t)  sin(OM_t)  0;
                    -sin(OM_t)  cos(OM_t)  0;
                        0          0       1];

            R1_i  = [ 1          0          0;
                    0       cos(i_t)   sin(i_t);
                    0      -sin(i_t)   cos(i_t)];

            RT_ELIO_to_Plane = R1_i * R3_OM;

            r1_plane = RT_ELIO_to_Plane * r1_current;
            r2_plane = RT_ELIO_to_Plane * r2_current;

            alpha1 = atan2(r1_plane(2), r1_plane(1));
            alpha2 = atan2(r2_plane(2), r2_plane(1));

            RT_Plane_to_ELIO = RT_ELIO_to_Plane';

            for om_tn_val = om_tn

                th_i_tn = alpha1 - om_tn_val;
                th_f_tn = alpha2 - om_tn_val;

                num_e = r2_n_current - r1_n_current;
                den_e = r1_n_current*cos(th_i_tn) - r2_n_current*cos(th_f_tn);
                e_t = num_e / den_e;

                if e_t < 0 || e_t >= 1, continue; end % Salto ellissi non valide

                num_a = r1_n_current * (1 + e_t*cos(th_i_tn));
                den_a = 1 - e_t^2;
                a_t = num_a / den_a;           
                rp_t = a_t * (1 - e_t);

                if rp_t <= R_sun + 1000, continue; end

                % Calcolo velocità (Logica par2car integrata ed efficiente)
                p_t = a_t * (1 - e_t^2);
                v_coeff = sqrt(mu/p_t);

                % Velocità nel frame perifocale (PF)
                v1t_pf = v_coeff * [-sin(th_i_tn); e_t + cos(th_i_tn); 0];
                v2t_pf = v_coeff * [-sin(th_f_tn); e_t + cos(th_f_tn); 0];

                R3_om = [ cos(om_tn_val)  sin(om_tn_val)  0;
                        -sin(om_tn_val)  cos(om_tn_val)  0;
                                0               0         1];

                R3_om_T = R3_om';

                v1t = RT_Plane_to_ELIO * (R3_om_T * v1t_pf);
                v2t = RT_Plane_to_ELIO * (R3_om_T * v2t_pf);     

                DeltaV = norm(v1t - v_i_current) + norm(v_f_current - v2t);
                    
                if DeltaV < best_dv

                    best_dv = DeltaV;
                    best_th_earth = th_i(j); 
                    best_th_nyx = th_f(k);
                    best_om_tn = om_tn_val;
                end
            end
        end
    end
end