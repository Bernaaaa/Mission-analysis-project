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


%% 1. CALCOLO PARAMETRI ORBITA DI PARTENZA E VINF

% Dati orbita di parcheggio alla fine dello scenario 1;
r_xi = -7090.590200;
r_yi = -5612.557300;
r_zi = 3948.902900;
v_xi = 5.698000;
v_yi = -5.995000;
v_zi = 1.710000;

% mu Terra
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

% Versore di V_inf in ECI (direzione asintotica dell'iperbole di uscita)
r_inf = v_inf_earth_eci / norm(v_inf_earth_eci);

% Imponendo la conservazione dell'energia meccanica, si può calcolare il semiasse maggiore dell'orbita iperbolica di arrivo:
ah = - mu / (norm(v_inf_earth_eci)^2);


%% TENTATIVO 1 - PROVO PERICENTRO ORBITALE DI PARTENZA COME PUNTO DI INIEZIONE
% Scelgo il pericentro dell'orbita di partenza come punto di iniezione, quindi rh è il raggio del pericentro dell'orbita di partenza e il raggio del punto di intersezione sull'ellittica
[rh_vec, v1] = par2car(a_i, e_i, i_i, OM_i, om_i, 0, mu);
% Prendo il versore di rh
rh_ver = rh_vec / norm(rh_vec);


%differenza angolare il raggio d'intersezione e la direzione asintotica
alpha = acos(rh_ver' * r_inf);


%% 2. CALCOLO PARAMETRI ORBITA IPERBOLICA DI ARRIVO E CONDIZIONE DI NON COPLANARITA'


% Raggio medio terrestre in km
r_earth = 6378.1363;

% torno a parametri kepleriani per determinare l'orbita
rh = norm(rh_vec);


x0 = [deg2rad(10); deg2rad(20); 2]; % guess iniziale scelta a posteriori guardando i risultati delle prime iterazioni

    %RISOLVO IL SISTEMA NON LINEARE DATO UN RH CON IL METODO DI NEWTON

    %imposto funzione
    f = @(x) [ ( ah * (1 - x(3)^2) ) / (1 + x(3) .* cos(x(2)))  - rh; ...
                                    cos(x(1)) + 1 / x(3); ...
                                    -x(1) + alpha + x(2)];

    % imposto metodo di newton
    fsol = fsolve(f, x0, optimoptions('fsolve', 'Display', 'final-detailed', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10 ));

    % ottengo parametri orbita iperbolica di uscita e li stampo a schermo
    theta_inf = fsol(1);
    theta_h = fsol(2);
    e_h = fsol(3);

    rp_h = ah * (1 - e_h^2);

    %check sulle condizioni che garantiscono un senso fisico alla soluzione (vd. pdf SC3_avanzato)
    if not (theta_inf > pi / 2 && theta_inf < pi  && rp_h > r_earth)
        warning('La soluzione trovata non soddisfa le condizioni di non coplanarità e di pericentro maggiore del raggio terrestre. Si consiglia di modificare la stima iniziale x0 e riprovare.');

    end


%stampo a schermo i risultati
fprintf('\n============================================================\n');
fprintf('| %-25s | %-15s |\n', 'Parametro', 'Valore');
fprintf('============================================================\n');

fprintf('| %-25s | %15.3f |\n', 'theta_inf [deg]', rad2deg(theta_inf));
fprintf('| %-25s | %15.3f |\n', 'theta_h [deg]',   rad2deg(theta_h));
fprintf('| %-25s | %15.6f |\n', 'e_h [-]',         e_h);
fprintf('| %-25s | %15.3f |\n', 'rh [km]',         rh);
fprintf('| %-25s | %15.3f |\n', 'v_inf [km/s]',    norm(v_inf_earth_eci));
fprintf('| %-25s | %15.3f |\n', 'a_h [km]',        ah);

fprintf('============================================================\n');


%% CARATTERIZZO COMPLETAMENTE IL PIANO DELL'ORBITA IPERBOLICA DI PARTENZA



% vettore momento angolare specifico dell'iperbole
%N.B potenzialmente h può avere qualunque verso : sia h che -h sono soluzioni accettabili e perpendicolare a r_inf, rh_ver (bisogna approfondire meglio questo punto)
h_ver = cross(rh_ver, r_inf) / norm(cross(rh_ver, r_inf)); 
h_H = sqrt(mu * ah * (1 - e_h^2)) .* h_ver; %vettore momento angolare specifico dell'iperbole


% vettore nodo ascendente dell'iperbole
k = [0; 0; 1];
N_H =  cross(k, h_H) / norm(cross(k, h_H)); 
N_H = N_H / norm(N_H); %dovrebbe essere già normalizzato, ma lo normalizzo comunque per sicurezza (better safe than sorry)

% inclinazione dell'iperbole
i_H = acos(h_H(3) / norm(h_H));


%ascensione retta del nodo ascendente
    if N_H(2) >= 0
        OM_H = acos(N_H(1));
    else
        OM_H = 2*pi - acos(N_H(1));
    end

% argomento del pericentro

% scelgo di usare il metodo definendo beta perché più veloce da implementare, ma si potrebbe tranquillmanente risolvere un altro sistema lineare (vd pdfSC3_avanzato)
cos_beta = dot(N_H, rh_ver);
sin_beta = dot(cross(N_H, rh_ver), (h_H / norm(h_H)) );

beta = atan2(sin_beta, cos_beta);
om_H = beta - theta_h;

% Calcolo vh
[rh1, vh] = par2car(ah, e_h, i_H, OM_H, om_H, theta_h, mu);


DV5 = norm(vh - v1);

% LA DV2 si calcola anologamente a quanto fatto nel caso coplanare
%% LAVORO DI OTTIMIZZAZIONE DELLA MANOVRA DI INIEZIONE

%chiaramente abbiamo sfiga e il pericentro fa cagare, quindi dobbiamo ottimizzare la manovra di iniezione per minimizzare la DV.
% ho fatto girare tutti i punti del'orbita di partenza per provarli uno a uno come punti di iniezione e usando un semplice grid search ho trovato il minimo

theta_low = 0;
theta_high = 2 * pi;
delta_theta = deg2rad(1); % passo di 1 grado

theta_val = theta_low:delta_theta:theta_high;
DV_values = zeros(size(theta_val));

for j = 1:length(theta_val)
    theta_i = theta_val(j);

    [rh_vec, v1] = par2car(a_i, e_i, i_i, OM_i, om_i, theta_i, mu);
    rh_ver = rh_vec / norm(rh_vec);

    rh = norm(rh_vec);

    alpha = acos(rh_ver' * r_inf);

    f = @(x) [ ( ah * (1 - x(3)^2) ) / (1 + x(3) .* cos(x(2)))  - rh; ...
                                    cos(x(1)) + 1 / x(3); ...
                                    -x(1) + alpha + x(2)];

    fsol = fsolve(f, x0, optimoptions('fsolve', 'Display', 'final-detailed', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10 ));

    theta_inf = fsol(1);
    theta_h = fsol(2);
    e_h = fsol(3);

    rp_h = ah * (1 - e_h^2);

    if theta_inf > pi / 2 && theta_inf < pi  && rp_h > r_earth
       
        % vettore momento angolare specifico dell'iperbole
        h_ver = cross(rh_ver, r_inf) / norm(cross(rh_ver, r_inf)); 
        h_H = sqrt(mu * ah * (1 - e_h^2)) .* h_ver; %vettore momento angolare specifico dell'iperbole
        
        control = true;

        while  control == true


        % vettore nodo ascendente dell'iperbole
        k = [0; 0; 1];
        N_H =  cross(k, h_H) / norm(cross(k, h_H)); 
        N_H = N_H / norm(N_H); %dovrebbe essere già normalizzato, ma lo normalizzo comunque per sicurezza (better safe than sorry)

        % inclinazione dell'iperbole
        i_H = acos(h_H(3) / norm(h_H));



        %ascensione retta del nodo ascendente
        if N_H(2) >= 0
            OM_H = acos(N_H(1));
        else
            OM_H = 2*pi - acos(N_H(1));
        end

         % argomento del pericentro

        % METODO 1 -  definendo beta
        cos_beta = dot(N_H, rh_ver);
        sin_beta = dot(cross(N_H, rh_ver), (h_H / norm(h_H)) );

        beta = atan2(sin_beta, cos_beta);
        om_H = beta - theta_h;

        % Calcolo vh
        [rh1, vh] = par2car(ah, e_h, i_H, OM_H, om_H, theta_h, mu);
        DV1 = norm(vh - v1);

        DV_values(j) = DV1;

        if dot(vh, v_inf_earth_eci) < 0
            h_H = - h_H; % se la direzione di vh è opposta a quella di v_inf vuol dire che abbiamo scelto il verso di percorrenza sbagliato dell'iperbole
        else
            control = false;
        end

    end
    else
        DV_values(j) = inf; % se la soluzione non soddisfa le condizioni di non coplanarità e di pericentro maggiore del raggio terrestre, assegno un valore infinito alla DV per escludere questa soluzione dall'ottimizzazione
    end
end

% Trovo il theta che minimizza la DV
[~, min_index] = min(DV_values);
optimal_theta = theta_val(min_index);
fprintf('\nIl theta ottimale per minimizzare la DV è: %.2f gradi\n', rad2deg(optimal_theta));
fprintf('\nLa DV minima ottenuta è: %.3f km/s\n', DV_values(min_index));

figure(1)
scatter(rad2deg(theta_val), DV_values);
xlabel('Theta [deg]');
ylabel('Delta-V [km/s]');
title('Delta-V in funzione di Theta');
grid on;

figure(2)
scatter3(cos(theta_val), sin(theta_val), DV_values);

DV5


