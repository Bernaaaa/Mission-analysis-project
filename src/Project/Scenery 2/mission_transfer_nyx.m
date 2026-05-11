function DV = mission_transfer_nyx(x, a_i, e_i, i_i, OM_i, om_i, a_f, e_f, i_f, OM_f, om_f, mu, R_sun)
    
    % Assigning the optimization variables
    th_i = x(1);
    th_f = x(2);
    om_tn_val = x(3);
    
    [r1, v_i] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu);
    [r2, v_f] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu); 
    r1_n = norm(r1);
    r2_n = norm(r2);
                    
    p_vect_r1r2 = cross(r1, r2);
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

    r1_plane = RT_ELIO_to_Plane * r1;
    r2_plane = RT_ELIO_to_Plane * r2;

    alpha1 = atan2(r1_plane(2), r1_plane(1));
    alpha2 = atan2(r2_plane(2), r2_plane(1));

    % Ensure that alpha2 is always greater than alpha1 (forward transfer)
    while alpha2 < alpha1
        alpha2 = alpha2 + 2*pi;
    end

    RT_Plane_to_ELIO = RT_ELIO_to_Plane';


    th_i_tn = alpha1 - om_tn_val;
    th_f_tn = alpha2 - om_tn_val;

    num_e = r2_n - r1_n;
    den_e = r1_n*cos(th_i_tn) - r2_n*cos(th_f_tn);
    e_t = num_e / den_e;

    num_a = r1_n * (1 + e_t*cos(th_i_tn));
    den_a = 1 - e_t^2;
    a_t = num_a / den_a;            
    rp_t = a_t * (1 - e_t); 

    if e_t < 0 || e_t >= 1 || rp_t <= R_sun + 1000 || a_t <= 0
        DV = 1e12; 
        return; 
    end

    p_t = a_t * (1 - e_t^2);
    v_coeff = sqrt(mu/p_t);

    v1t_pf = v_coeff * [-sin(th_i_tn); e_t + cos(th_i_tn); 0];
    v2t_pf = v_coeff * [-sin(th_f_tn); e_t + cos(th_f_tn); 0];

    R3_om = [ cos(om_tn_val)  sin(om_tn_val)  0;
            -sin(om_tn_val)  cos(om_tn_val)  0;
                    0               0         1];

    R3_om_T = R3_om';

    v1t = RT_Plane_to_ELIO * (R3_om_T * v1t_pf);
    v2t = RT_Plane_to_ELIO * (R3_om_T * v2t_pf);     

    DV = norm(v1t - v_i) + norm(v_f - v2t);
     
end