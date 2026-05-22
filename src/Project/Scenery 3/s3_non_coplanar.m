clc
clear
close all

%% 0. PATH SETUP
% -------------------------------------------------------------------------
currentDir = fileparts(mfilename('fullpath'));
funcFolder = fullfile(currentDir, '..', '..', 'lab');

if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
else
    warning('Warning: Functions folder not found at: %s', funcFolder);
end



% Dati orbita di parcheggio alla fine dello scenario 1;
r_xi = -7090.590200;
r_yi = -5612.557300;
r_zi = 3948.902900;
v_xi = 5.698000;
v_yi = -5.995000;
v_zi = 1.710000;

mu = 398600.4418;


ri = [r_xi; r_yi; r_zi];
vi = [v_xi; v_yi; v_zi];


[a_i, e_i, i_i, OM_i, om_i, th_i] = car2par(ri, vi, mu);


%Calcolo Vinf con i risultati dello scenario 2;
%N.B. I dati di scenario 2 sono in eclittica, mentre quelli di scenario 1 sono in ECI, quindi è necessario convertire i dati di scenario 2 in ECI prima di calcolare Vinf.

V_Earth_Depart =  [  2.9766,  29.5101,  -0.0026];
V_Trans_Depart =  [  3.1273,  32.4459,   0.5897];

v_inf_earth = V_Trans_Depart - V_Earth_Depart;

% Rotazione da eclittica a ECI

epsilon = deg2rad(23.45);% [rad] Obliquity of the ecliptic

% Matrice di rotazione per convertire da eclittica a ECI
T = [1, 0, 0; 0, cos(epsilon), sin(epsilon); 0, - sin(epsilon), cos(epsilon)];


v_inf_earth_eci = T' * v_inf_earth';

% Versore di V_inf in ECI
r_inf = v_inf_earth_eci / norm(v_inf_earth_eci);

% Imponendo la conservazione dell'energia meccanica, si può calcolare il semiasse maggiore dell'orbita iperbolica di arrivo:
ah = - mu / (norm(v_inf_earth_eci)^2);

% Scelgo il pericentro dell'orbita di partenza come punto di iniezione, quindi rh è il raggio del pericentro dell'orbita di partenza, che è anche il raggio del pericentro dell'orbita iperbolica di arrivo, poiché i due orbiti sono tangenti in quel punto.
rh_ver = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);


% Prendo il versore di rh
rh_ver = rh_ver / norm(rh_ver);


alpha = acos(rh_ver' * r_inf);

control = true;


% Delta arbitrario per aumentare il raggio del pericentro dell'orbita iperbolica di arrivo, in modo da trovare un punto di iniezione che soddisfi la condizione di non coplanarità tra l'orbita di partenza e l'orbita iperbolica di arrivo.
delta = 10;

% Raggio medio terrestre in km
r_earth = 6378.1363;

% Inizializzo rh con un valore che è sicuramente inferiore al raggio del pericentro dell'orbita di partenza, in modo da garantire che il punto di iniezione sia all'interno dell'orbita di partenza.
rh = r_earth + 300;


x0 = [deg2rad(120); deg2rad(20); 2];
n_iter = 0;

while control == true

    n_iter = n_iter + 1;

    rh = rh + delta;

    %RISOLVO IL SISTEMA NON LINEARE DATO UN RH CON IL METODO DI NEWTON
    f = @(x) [ ( ah * (1 - x(3)^2) ) / (1 + x(3) .* cos(x(1)))  - rh; ...
                                    cos(x(1)) + 1 / x(3); ...
                                    x(1) - alpha - x(2)];

    fsol = fsolve(f, x0, optimoptions('fsolve', 'Display', 'off'));
    theta_inf = fsol(1);
    theta_h = fsol(2);
    e_h = fsol(3);

    if theta_inf < pi && theta_inf > pi / 2
        control = false;

    end

fprintf('\npensando2');



end


fprintf('\n il metodo converge dopo %d iterazioni\n', n_iter)

fprintf('\n============================================================\n');
fprintf('| %-25s | %-15s |\n', 'Parametro', 'Valore');
fprintf('============================================================\n');

fprintf('| %-25s | %15.3f |\n', 'theta_inf [deg]', rad2deg(theta_inf));
fprintf('| %-25s | %15.3f |\n', 'theta_h [deg]',   rad2deg(theta_h));
fprintf('| %-25s | %15.6f |\n', 'e_h [-]',         e_h);
fprintf('| %-25s | %15.3f |\n', 'rp [km]',         rh);
fprintf('| %-25s | %15.3f |\n', 'v_inf [km/s]',    norm(v_inf_earth_eci));
fprintf('| %-25s | %15.3f |\n', 'a_h [km]',        ah);

fprintf('============================================================\n');


%% CARATTERIZZO COMPLETAMENTE IL PIANO DELL'ORBITA IPERBOLICA DI PARTENZA

% vettore momento angolare specifico dell'iperbole
h = cross(rh_ver .* rh, v_inf_earth_eci);

%direzione del vettore momento angolare
h_ver = h / norm(h);

% vettore nodo ascendente dell'iperbole
k = [0; 0; 1];
N = cross(k, h / norm(h));
N = N / norm(N); %dovrebbe essere già normalizzato, ma lo normalizzo comunque per sicurezza (better safe than sorry)

% inclinazione dell'iperbole
i = acos(h(3));

 %ascensione retta del nodo ascendente
    if N(2) >= 0
        OM = acos(N(1));
    else
        OM = 2*pi - acos(N(1));
    end

% argomento del pericentro

% METODO 1 -  definendo beta
cos_beta = dot(N, rh_ver);
sin_beta = dot(cross(N, rh_ver), h / norm(h));

omH = atan2(sin_beta, cos_beta);


[~, vh] = par2car(ah, e_h, i, OM, omH, theta_h, mu);
[~, v1] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);   

DV1 = norm(vh - v1);

% LA DV2 si calcola anologamente a quanto fatto nel caso coplanare


        

