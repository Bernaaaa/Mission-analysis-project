%Strategy 2: Bielliptic Transfer + Plane Change at Apogee + Clean-up at Final Orbit
% The objective is to maximize the efficiency of the transfer by performing the plane change at the apogee of a 
% highly elliptical transfer orbit, where the velocity is lower, and then doing a final clean-up maneuver 
% to adjust the argument of pericenter.
close all; clear; clc;

%% --- SETUP PERCORSI UNIVERSALE ---
currentDir = fileparts(mfilename('fullpath'));

funcFolder = fullfile(currentDir, '..', '..', 'lab' );

if exist(funcFolder, 'dir')
    addpath(genpath(funcFolder));
else
    warning('Attenzione: La cartella delle funzioni non è stata trovata in: %s', funcFolder);
end

%% --- CALCOLO GRIGLIA DI VALORI ---

mu = 398600.4418;

a = 24400.00;
e = 0.728300;
i = 0.104700;
OM = 2.361000;
om = 3.107000;
th = 2.135000;

r_xf = -7090.590200;
r_yf = -5612.557300;
r_zf = 3948.902900;
v_xf = 5.698000;
v_yf = -5.995000;
v_zf = 1.710000;

rf = [r_xf; r_yf; r_zf];
vf = [v_xf; v_yf; v_zf];

rfm = norm(rf,2);
vfm = norm(vf,2);

[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rf, vf, mu);

ra_range = linspace(50000, 315000, 500);
 % Apocenter of the transfer orbit (chosen to be very high for a more efficient bielliptic transfer)
e_range = linspace(0.01, 0.95, 500); 

DV_val = zeros(length(ra_range), length(e_range));
DV_val1 = DV_val;
DV_val2 = DV_val;
DV_val3 = DV_val;
DV_val4 = DV_val;
DV_val5 = DV_val;
DV_val6 = DV_val;

best_dv = inf;
best_ra = 0;
best_e = 0;


idx_ra = 0;
idx_e = 0;

dist_ra = ra_range(2) - ra_range(1);
ra_range = [ra_range, 315000 + dist_ra : dist_ra : 350000 + dist_ra]; % Estensione del range di apocentri per includere anche valori più alti


for r_test = ra_range
    idx_ra = idx_ra + 1;
    idx_e = 0;

    for e_test = e_range

        idx_e = idx_e + 1;
        
        a_test = r_test / (1 + e_test);
        
        rp_test = a_test * (1 - e_test);
        if rp_test < 6378 + 200 % Impossibile: pericentro troppo basso (considerando un'altitudine minima di 200 km)
            DV_val(idx_ra, idx_e) = NaN; % Impossibile: pericentro troppo basso
            continue;
        end
        
        try 
            [DV1a, DV1b] = bitangentTransfer(a, e, a_test, e_test, 'pa', mu);
            [DVP, om_p] = changeOrbitalPlaneNoPrint(a_test, e_test, i, OM, om, i_f, OM_f, mu); % CAMBIA E!!!
            [DVom] = changePericenterArg(a_test, e_test, om_p, om_f, mu);
            [DV2a, DV2b] = bitangentTransfer(a_test, e_test, a_f, e_f, 'ap', mu);
            
            DV_tot = abs(DV1a) + abs(DV1b) + abs(DVP) + abs(DVom) + abs(DV2a) + abs(DV2b);
   

            DV_val(idx_ra, idx_e) = DV_tot;
            DV_val1(idx_ra, idx_e) = abs(DV1a);
            DV_val2(idx_ra, idx_e) = abs(DV1b);
            DV_val3(idx_ra, idx_e) = abs(DVP);
            DV_val4(idx_ra, idx_e) = abs(DVom);
            DV_val5(idx_ra, idx_e) = abs(DV2a);
            DV_val6(idx_ra, idx_e) = abs(DV2b);
            
            % --- Salvataggio del Record ---
            if DV_tot < best_dv && r_test <= 315000
                best_dv = DV_tot;
                best_ra = r_test;
                best_e = e_test;
            end
        catch
            continue; % Se la geometria è impossibile, passa alla prossima
        end
    end
end

%% PLOT SCENARIO 1: DELTA-V MAP

L = 500; %indice t.c. ra = 315000km



%colore linea principale
colore_main = [0.10, 0.10, 0.10]; 

%colore contributi singoli
colori_contributi = [
    0.45, 0.62, 0.81;  
    1.00, 0.62, 0.29;  
    0.40, 0.75, 0.50;  
    0.93, 0.40, 0.36;  
    0.68, 0.55, 0.79;  
    0.93, 0.59, 0.79   
];


%% GRAFICO 1: DELTA V TO RA

figure()

%il bordo del dominio è posto all'inizio in modo da trovarsi sullo sfondo rispetto alle altre linee
p9 = plot([best_ra ./ 1e+5, best_ra ./ 1e+5], [0, 4.5], '--', 'Color', [0.35, 0.40, 0.45], ...
    'DisplayName', 'Upper Bound', 'LineWidth', 2); 
hold on

%DV totale
p1 = plot(ra_range ./ 1e+5, DV_val(:, 1), 'Color', colore_main, 'LineWidth', 2, 'DisplayName', '$\Delta V$ total');

%Singole manovre
p2 = plot(ra_range ./ 1e+5, DV_val1(:, 1), 'Color', colori_contributi(1,:), 'LineWidth', 2, 'DisplayName', 'Departing Bitangent 1° impulse');
p3 = plot(ra_range ./ 1e+5, DV_val2(:, 1), 'Color', colori_contributi(2,:), 'LineWidth', 2, 'DisplayName', 'Departing Bitangent 2° impulse');
p4 = plot(ra_range ./ 1e+5, DV_val3(:, 1), 'Color', colori_contributi(3,:), 'LineWidth', 2, 'DisplayName', 'Plane Change');
p5 = plot(ra_range ./ 1e+5, DV_val4(:, 1), 'Color', colori_contributi(4,:), 'LineWidth', 2, 'DisplayName', 'Argument Pericenter Adjustment');
p6 = plot(ra_range ./ 1e+5, DV_val5(:, 1), 'Color', colori_contributi(5,:), 'LineWidth', 2, 'DisplayName', 'Arriving Bitangent 3° impulse');
p7 = plot(ra_range ./ 1e+5, DV_val6(:, 1), 'Color', colori_contributi(6,:), 'LineWidth', 2, 'DisplayName', 'Arriving Bitangent 4° impulse');

%punto ottimale
p8 = plot(best_ra ./ 1e+5, best_dv, ...
    'LineStyle', 'none', ...
    'Marker', 'o', ...             
    'MarkerSize', 9, ...           
    'LineWidth', 2.0, ...           % Spessore del bordo
    'MarkerEdgeColor', [0.2 0.2 0.2], ... 
    'MarkerFaceColor', [1.000 0.600 0.000], ... 
    'DisplayName', 'Optimal Point');

xlabel('$r_a$  $[10^5 km]$', 'FontSize', 12);
ylabel('$\Delta V [km/s]$', 'FontSize', 12);


lgd = legend([p1, p2, p3, p4, p5, p6, p7, p8, p9],...
    'Location', 'northeastoutside');
    lgd.FontSize = 10;
    lgd.Title.String = '$\Delta V$ contributions for each impulse';
    lgd.Title.Interpreter = 'latex';
    lgd.Box = 'on';
    
title('$\Delta V(e = 0.01, r_a)$', 'FontSize', 14);
grid on; axis tight;

%% GRAFICO 2: DELTA V TO E

figure()

plot(e_range, DV_val(L, :), 'Color', colore_main, 'LineWidth', 2, 'DisplayName', '$\Delta V$ total');
hold on

plot(e_range, DV_val1(L, :), 'Color', colori_contributi(1,:), 'LineWidth', 2, 'DisplayName', 'Departing Bitangent 1° impulse');
plot(e_range, DV_val2(L, :), 'Color', colori_contributi(2,:), 'LineWidth', 2, 'DisplayName', 'Departing Bitangent 2° impulse');
plot(e_range, DV_val3(L, :), 'Color', colori_contributi(3,:), 'LineWidth', 2, 'DisplayName', 'Plane Change');
plot(e_range, DV_val4(L, :), 'Color', colori_contributi(4,:), 'LineWidth', 2, 'DisplayName', 'Argument Pericenter Adjustment.');
plot(e_range, DV_val5(L, :), 'Color', colori_contributi(5,:), 'LineWidth', 2, 'DisplayName', 'Arriving Bitangent 1° impulse');
plot(e_range, DV_val6(L, :), 'Color', colori_contributi(6,:), 'LineWidth', 2, 'DisplayName', 'Arriving Bitangent 2° impulse');

plot(best_e, best_dv, ...
    'LineStyle', 'none', ...
    'Marker', 'o', ...              
    'MarkerSize', 9, ...          
    'LineWidth', 2.0, ...           
    'MarkerEdgeColor', [0.2 0.2 0.2], ... 
    'MarkerFaceColor', [1.000 0.600 0.000], ... 
    'DisplayName', 'Optimal Point');

xlabel('$e$', 'FontSize', 12);
ylabel('$\Delta V [km/s]$', 'FontSize', 12);

lgd = legend('show', 'Location', 'northeastoutside');
    lgd.FontSize = 10;
    lgd.Title.String = '$\Delta V$ contributions for each impulse';
    lgd.Title.Interpreter = 'latex';
    lgd.Box = 'on';

title('$\Delta V(e, r_a = 315000km)$', 'FontSize', 14);
grid on; axis tight;
