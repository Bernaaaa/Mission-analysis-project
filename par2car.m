function [rr, vv] = par2car(a, e, i, OM, om, th, mu)
    %La funzione par2car calcola il vettore di stato [r, v] in coordinate ECI a partire dai parametri orbitali
    %
    % [rr, vv] = par2car(a, e, i, OM, om, th, mu)
    %
    %INPUT:
    % -----------------------------------------------------------------
    % a - semiasse maggiore in [km]
    % e - eccentricità
    % i - inclinazione in [rad]
    % OM - RAAN in [rad]
    % om - argomento del pericentro in [rad]
    % th - anomalia vera in [rad]
    % mu - parametro gravitazionale standard in [km^3/s^2]

    % -----------------------------------------------------------------
    %OUTPUT:
    % rr - vettore posizione [x, y, z] in [km]
    % vv - vettore velocità [vx, vy, vz] in [km/s]

    p = a * (1 - e^2); %parametro orbitale

    r = p / (1 + e * cos(th)); %raggio

    %posizione nel piano orbitale
    x_orb = r * cos(th);
    y_orb = r * sin(th);
    z_orb = 0;

    rr_pf = [x_orb; y_orb; z_orb];

    %velocità nel piano orbitale
    vx_orb = -sqrt(mu/p) * sin(th);
    vy_orb = sqrt(mu/p) * (e + cos(th));
    vz_orb = 0;

    vv_pf = [vx_orb; vy_orb; vz_orb];

    %matrice di rotazione per passare dal piano orbitale al sistema di riferimento inerziale
    R3_OM = [cos(OM), sin(OM), 0;
            -sin(OM), cos(OM), 0; 
            0, 0, 1];
    R1_i = [1, 0, 0; 
            0, cos(i), sin(i); 
            0, -sin(i), cos(i)];

    R3_om = [cos(om), sin(om), 0; 
            -sin(om), cos(om), 0; 
            0, 0, 1];

    T = R3_OM' * R1_i' * R3_om';

    rr = T * rr_pf;
    vv = T * vv_pf;
end