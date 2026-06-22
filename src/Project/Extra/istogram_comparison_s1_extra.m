%% ISTOGRAMMA DELLE STRATEGIE DI TRASFERIMENTO
clear; clc; close all;
%% DATI STRATEGIA 1: STANDARD

%eseguo scenery_1_standard per ottenere i dati della strategia 1
scenery1_standard; 
clc
close all

DV_pl_standard = abs( DeltaVP_BLT);
DV_BT_standard = abs(DeltaV1BT) + abs(DeltaV2BT); 
DV_cp_standard = abs(DeltaVCP);
DV_total_standard = abs(DeltaV_total);
DT_standard = DeltaT_TOT;

% salvo su file temporaneo (per evitare che vengano cancellati dei clear in testa agli altri file)
save('temp_strat1.mat', 'DV_pl_standard', 'DV_BT_standard', 'DV_cp_standard', 'DV_total_standard', 'DT_standard');
clear

%% DATI STRATEGIA 2: DELTA T
% Anche se questo script fa un "clear", noi abbiamo i dati al sicuro sul disco
scenery1_deltaT; 
clc
close all

DV_pl_Delta_T  = abs(dv_plane);
DV_BT_Delta_T  = abs(DeltaV1BT) + abs(DeltaV2BT); 
DV_cp_Delta_T  = abs(DeltaVCP); 
DV_total_Delta_T = abs(DeltaV_total);
DT_Delta_T = DeltaT_TOT; 

save('temp_strat2.mat', 'DV_pl_Delta_T', 'DV_BT_Delta_T', 'DV_cp_Delta_T', 'DV_total_Delta_T', 'DT_Delta_T');
clear

%% DATI STRATEGIA 3: DELTA V
scenery1_deltaV;
clc
close all

DV_pl_deltaV = abs (DeltaVP_BLT);
DV_cp_deltaV = abs (DeltaV_omega_blt);
DV_BT_deltaV = (abs(DeltaV1BT) + abs(DeltaV2BT)) + (abs(DeltaV3BT) + abs(DeltaV4BT));
DV_total_deltaV = abs (DeltaV_Totale_BLT);
DT_DeltaV = DeltaT_TOT;

save('temp_strat3.mat', 'DV_pl_deltaV', 'DV_cp_deltaV', 'DV_BT_deltaV', 'DV_total_deltaV', 'DT_DeltaV');

clear

%% COSTRUZIONE MATRICE E CANCELLAZIONE FILE TEMPORANEI

% Carico i dati salvati nei file temporanei
load('temp_strat1.mat');
load('temp_strat2.mat');
load('temp_strat3.mat')

% Creo le colonne della matrice dei dati: ognuna ordinate con le sue manovre per ogni strategia
colonna_Standard = [DV_BT_standard(1); DV_pl_standard(1); DV_cp_standard(1); DV_total_standard(1)];
colonna_DeltaT   = [DV_BT_Delta_T(1); DV_pl_Delta_T(1); DV_cp_Delta_T(1); DV_total_Delta_T(1)];
colonna_DeltaV   = [DV_BT_deltaV(1); DV_pl_deltaV(1); DV_cp_deltaV(1); DV_total_deltaV(1)];

data = [colonna_Standard, colonna_DeltaT, colonna_DeltaV];

% elimino i temporanei
delete('temp_strat1.mat');
delete('temp_strat2.mat');
delete('temp_strat3.mat');

%% CREAZIONE GRAFICO

%imposto la figura con sfondo bianco
figure('Color', [1 1 1]); 
hold on; 


% Creo le barre dell'istogramma
hBar = bar(data, 'grouped', 'EdgeColor', 'none'); 

% Selezione della palette di colori:
hBar(1).FaceColor = [0.85, 0.45, 0.25]; % Standard Strategy (Blu)
hBar(2).FaceColor = [0.20, 0.45, 0.65]; % Delta T Strategy (Arancione)
hBar(3).FaceColor = [0.25, 0.60, 0.45]; % Delta V Strategy (Verde)


%configuro l'asse x
set(gca, 'XTick', 1:4); %imposto il numero di raggruppamenti di colonne

set(gca, 'XTickLabel', {'$\Delta V_{BT}$', '$\Delta V_{plane}$', '$\Delta V_{CP}$', '$\Delta V_{Total}$'}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11); %nomino l'asse

ylabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('Comparison of Transfer Strategies', 'FontSize', 13, 'FontWeight', 'bold');

%% 8. LEGENDA (Legata solo alle 3 colonne)
% Passando explicitamente 'hBar' (che ha lunghezza 3), la legenda mostrerà solo 3 elementi
lgd = legend(hBar, {'Standard Strategy', '$\Delta T$ Strategy', '$\Delta V$ Strategy'}, ...
    'Location', 'eastoutside', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.String = 'Transfer Strategies';
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';

hold off;


%% CREAZIONE GRAFICO ORIZZONTALE
figure('Color', [1 1 1]); 
hold on; 
grid on;

blu      = [31, 119, 180] / 255;
arancio  = [255, 127, 14] / 255;
verde    = [44, 160, 44]  / 255;

colonna_Standard = [DV_total_standard(1); DV_cp_standard(1); DV_pl_standard(1); DV_BT_standard(1)];
colonna_DeltaT = [DV_total_Delta_T(1); DV_cp_Delta_T(1); DV_pl_Delta_T(1); DV_BT_Delta_T(1)]
colonna_DeltaV = [DV_total_deltaV(1); DV_cp_deltaV(1); DV_pl_deltaV(1);DV_BT_deltaV(1);]

data = [colonna_Standard, colonna_DeltaT, colonna_DeltaV];

% --- NOTA IL COMANDO 'barh' AL POSTO DI 'bar' ---
hBar = barh(data, 'grouped', 'EdgeColor', 'none'); 

% Applichiamo i 3 colori alle 3 colonne (Strategie)
hBar(1).FaceColor = arancio; % Strategia 1 (arancione)
hBar(2).FaceColor = blu; % Strategia 2 (blu)
hBar(3).FaceColor = verde; %Srategia 3 (Verde)

%% 6. CONFIGURAZIONE ASSI INVERTITA
% Ora le etichette vanno sull'asse Y (YTick e YTickLabel)
set(gca, 'YTick', 1:4);
set(gca, 'YTickLabel', {'$\Delta V_{TOT}$', '$\Delta V_{per}$', '$\Delta V_{plane}$', '$\Delta V_{BT}$'}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Invertiamo anche i titoli degli assi X e Y rispetto a prima
xlabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('$\Delta V$ Strategy Comparison', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

%% 7. LEGENDA (Rimane a destra, agganciata correttamente)
lgd = legend(hBar, {'First Strategy', 'Second Strategy', 'Third Strategy'}, ...
    'Location', 'northeastoutside', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.String = 'Transfer Strategies';
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';


hold off;

%% ISTOGRAMMA COMPARISON DELTA V & DELTA T FOR EACH STRATEGY



DV_data = [DV_total_standard(1); DV_total_Delta_T(1); DV_total_deltaV(1)];
DT_data = [DT_standard(1); DT_Delta_T(1); DT_DeltaV(1)];

data = [DV_data, DT_data];
data(:, 2) = data(:, 2) ./ 10^5;
figure('Color', [1 1 1]); 
hold on; 
grid on;
ax = gca; 
ax.GridColor = [0.9 0.9 0.9];

% Disegniamo le barre raggruppate
hBar = bar(data, 'grouped', 'EdgeColor', 'none'); 

hBar(1).FaceColor = [0.0, 0.4, 1.0]; % Colonna 1: Delta V (Blu)
hBar(2).FaceColor = [1.0, 0.1, 0.1]; % Colonna 2: Delta T (Rosso)

% Configurazione asse X
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Strategy 1', 'Strategy 2', 'Strategy 3'}, ...
         'FontSize', 11); % Rimosso il $ del latex se usi testo normale, per una resa più pulita

% --- MODIFICA ETICHETTA ASSE Y ---
% Ora l'asse Y spiega entrambe le unità di misura convertite
ylabel('$\Delta V$ [km/s] \quad / \quad $\Delta T$ [$10^5$s]', 'Interpreter', 'latex', 'FontSize', 12);


%% 8. LEGENDA (Sincronizzata con le colonne)
lgd = legend(hBar, {'$\Delta V$', '$\Delta T$'}, ...
    'Location', 'northwest', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';

hold off;