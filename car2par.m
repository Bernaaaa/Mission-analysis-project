function [a, e, i, OM, om, th] = car2par(rr, vv, mu)
    % Dato un vettore di stato [r, v] in coordinate ECI, calcola i parametri orbitali

    % [a, e, i, OM, om, th] = car2par(rr, vv, mu)

    %---------------------------------------------------------------
    %INPUT:

    % rr - vettore posizione [x, y, z] in km
    % vv - vettore velocità [vx, vy, vz] in km/s
    % mu - parametro gravitazionale standard in km^3/s^2

    %---------------------------------------------------------------
    %OUTPUT:
    % a - semiasse maggiore in km
    % e - eccentricità
    % i - inclinazione in rad
    % omega - argomento del pericentro in rad
    % w - RAAN in rad
    % theta - anomalia vera in rad

    r = norm(rr);
    v = norm(vv);

    %Semiasse maggiore
    a = 1 / (2/r - v^2/mu); 

    %momento angolare
    hh = cross(rr, vv);

    %eccentricità
    ee =  cross(vv, hh)/mu - rr/r;
    e = norm(ee);

    %inclinazione
    i = acos(hh(3)/norm(hh));

    %versore k in direzione z
    k = [0; 0; 1];

    %Linea dei nodi
    N = cross(k, hh) / norm(cross(k, hh)); 

    %ascensione retta del nodo ascendente
    if N(2) >= 0
        OM = acos(N(1));
    else
        OM = 2*pi - acos(N(1));
    end

    %anomalia al pericentro
    if ee(3) >= 0
        om = acos(dot(N, ee)/e);
    else
        om = 2*pi - acos(dot(N, ee)/e);
    end

    %velocità radiale
    vr = dot(rr, vv) / r;

    if vr >= 0 %se la velocità radiale è positiva, allora l'anomalia vera è tra 0 e pi
        th = acos(dot(ee, rr)/(e*r));
    else
        th = 2*pi - acos(dot(ee, rr)/(e*r));
    end





    