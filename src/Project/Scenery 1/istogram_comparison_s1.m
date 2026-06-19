%% ISTOGRAMMA DELLE STRATEGIE DI TRASFERIMENTO
clear; clc; close all;

%% 1. STRATEGIA 1: STANDARD
scenery1_standard; 
clc
close all

% Calcolo e salvataggio immediato
DV_pl_standard = abs( DeltaVP_BLT);
DV_BT_standard = abs(DeltaV1BT) + abs(DeltaV2BT); 
DV_cp_standard = abs(DeltaVCP);
DV_total_standard = abs(DeltaV_total);
DT_standard = DeltaT_TOT 

% SALVAGGIO DI SICUREZZA SU DISCO
save('temp_strat1.mat', 'DV_pl_standard', 'DV_BT_standard', 'DV_cp_standard', 'DV_total_standard', 'DT_standard');


clear


%% 2. STRATEGIA 2: DELTA T
% Anche se questo script fa un "clear", noi abbiamo i dati al sicuro sul disco
scenery1_deltaT; 
clc
close all

DV_pl_Delta_T  = abs(dv_plane);
DV_BT_Delta_T  = abs(DeltaV1BT) + abs(DeltaV2BT); 
DV_cp_Delta_T  = abs(DeltaVCP); 
DV_total_Delta_T = abs(DeltaV_total);
DT_Delta_T = DeltaT_TOT 

% SALVAGGIO DI SICUREZZA SU DISCO
save('temp_strat2.mat', 'DV_pl_Delta_T', 'DV_BT_Delta_T', 'DV_cp_Delta_T', 'DV_total_Delta_T', 'DT_Delta_T');


clear

%% 3. STRATEGIA 3: DELTA V
scenery1_deltaV;
clc
close all

DV_pl_deltaV = abs (DeltaVP_BLT);
DV_cp_deltaV = abs (DeltaV_omega_blt);
DV_BT_deltaV = (abs(DeltaV1BT) + abs(DeltaV2BT)) + (abs(DeltaV3BT) + abs(DeltaV4BT));
DV_total_deltaV = abs (DeltaV_Totale_BLT);
DT_DeltaV = DeltaT_TOT 

save('temp_strat3.mat', 'DV_pl_deltaV', 'DV_cp_deltaV', 'DV_BT_deltaV', 'DV_total_deltaV', 'DT_DeltaV');

clear


%% 4. IL MOMENTO DEL RECUPERO DATI
% Puliamo tutto il caos generato dagli script e ricarichiamo solo il necessario


load('temp_strat1.mat');
load('temp_strat2.mat');
load('temp_strat3.mat')

%% 5. COSTRUZIONE MATRICE E CANCELLAZIONE FILE TEMPORANEI
% Creiamo prima i vettori colonna per ogni strategia per essere sicuri al 100%
colonna_Standard = [DV_BT_standard(1); DV_pl_standard(1); DV_cp_standard(1); DV_total_standard(1)];
colonna_DeltaT   = [DV_BT_Delta_T(1); DV_pl_Delta_T(1); DV_cp_Delta_T(1); DV_total_Delta_T(1)];
colonna_DeltaV   = [DV_BT_deltaV(1); DV_pl_deltaV(1); DV_cp_deltaV(1); DV_total_deltaV(1)];

% Pulizia dei file temporanei creati per non lasciare sporco sul PC
delete('temp_strat1.mat');
delete('temp_strat2.mat');
delete('temp_strat3.mat');

% Assembliamo la matrice finale: 4 righe (contributi) e 3 colonne (strategie)
data = [colonna_Standard, colonna_DeltaT, colonna_DeltaV];

%% 6. CREAZIONE GRAFICO
figure('Color', [1 1 1]); 
hold on; 
grid on;
ax = gca; 
ax.GridColor = [0.9 0.9 0.9];

% Disegniamo le barre raggruppate
hBar = bar(data, 'grouped', 'EdgeColor', 'none'); 

% Applichiamo i 3 colori alle 3 colonne (Strategie)
hBar(1).FaceColor = [0.85, 0.45, 0.25]; % Standard Strategy (Blu)
hBar(2).FaceColor = [0.20, 0.45, 0.65]; % Delta T Strategy (Arancione)
hBar(3).FaceColor = [0.25, 0.60, 0.45]; % Delta V Strategy (Verde)

%% 7. CONFIGURAZIONE ASSI (4 Gruppi esatti)
% Forziamo l'asse X ad avere esattamente 4 tacche, una per ogni riga di 'data'
set(gca, 'XTick', 1:4);
set(gca, 'XTickLabel', {'$\Delta V_{BT}$', '$\Delta V_{plane}$', '$\Delta V_{CP}$', '$\Delta V_{Total}$'}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

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
ax = gca; 
ax.GridColor = [0.9 0.9 0.9];

% --- NOTA IL COMANDO 'barh' AL POSTO DI 'bar' ---
hBar = barh(data, 'grouped', 'EdgeColor', 'none'); 

% Applichiamo i 3 colori alle 3 colonne (Strategie)
hBar(1).FaceColor = [0.20, 0.45, 0.65]; % Standard Strategy (Blu)
hBar(2).FaceColor = [0.85, 0.45, 0.25]; % Delta T Strategy (Arancione)

%% 6. CONFIGURAZIONE ASSI INVERTITA
% Ora le etichette vanno sull'asse Y (YTick e YTickLabel)
set(gca, 'YTick', 1:4);
set(gca, 'YTickLabel', {'$\Delta V_{BT}$', '$\Delta V_{plane}$', '$\Delta V_{CP}$', '$\Delta V_{Total}$'}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Invertiamo anche i titoli degli assi X e Y rispetto a prima
xlabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('Comparison of Transfer Strategies', 'FontSize', 13, 'FontWeight', 'bold');

%% 7. LEGENDA (Rimane a destra, agganciata correttamente)
lgd = legend(hBar, {'Standard Strategy', '$\Delta T$ Strategy', '$\Delta V$ Strategy'}, ...
    'Location', 'eastoutside', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.String = 'Transfer Strategies';
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';

hold off;

%% ISTOGRAMMA COMPARISON DELTA V & DELTA T FOR EACH STRATEGY


figure('Color', [1 1 1]); 
hold on; 
grid on;

% Assicuriamoci che siano vettori colonna puliti
DV_data = [DV_total_standard(1); DV_total_Delta_T(1); DV_total_deltaV(1)];
DT_data = [DT_standard(1); DT_Delta_T(1); DT_DeltaV(1)];

% Coordinate base dell'asse X per le 3 strategie
x_coor = [1, 2, 3]; 

% --- PARAMETRI GEOMETRICI CRUCIALI ---
width_bar = 0.25;  % Spessore di ogni singola barra
offset    = 0.125;  % Distanza dal centro (totale separazione = 0.32)

%% --- ASSE Y SINISTRA: DELTA V (Barra Blu spostata a sinistra) ---
yyaxis left
hBar1 = bar(x_coor - offset, DV_data, width_bar, 'EdgeColor', 'none', 'FaceColor', [0.20, 0.45, 0.65]);
ylabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
ax = gca;
ax.YColor = [0, 0, 0]; % Numeri asse sinistro blu

%% --- ASSE Y DESTRA: DELTA T (Barra Arancione spostata a destra) ---
yyaxis right
hBar2 = bar(x_coor + offset, DT_data, width_bar, 'EdgeColor', 'none', 'FaceColor', [0.85, 0.45, 0.25]);
ylabel('$\Delta T$ [s]', 'Interpreter', 'latex', 'FontSize', 12);
ax.YColor = [0, 0, 0]; % Numeri asse destro arancione

%% --- CONFIGURAZIONE ASSI ---
ax.GridColor = [0.9 0.9 0.9];
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Strategy 1', 'Strategy 2', 'Strategy 3'}, 'FontSize', 11);

% Allarghiamo leggermente i limiti dell'asse X per non far toccare ai grafici i bordi del riquadro
xlim([0.5, 3.5]); 

title('Comparison of Transfer Strategies (Dual Axes)', 'FontSize', 13, 'FontWeight', 'bold');

%% --- LEGENDA ---
lgd = legend([hBar1, hBar2], {'$\Delta V$ (Left Axis)', '$\Delta T$ (Right Axis)'}, ...
    'Location', 'eastoutside', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Box = 'on';

hold off;

%% ALTERNATIVA CON DIVERSE SCALE


data = [DV_data, DT_data];
data(:, 2) = data(:, 2) ./ 10^5;
figure('Color', [1 1 1]); 
hold on; 
grid on;
ax = gca; 
ax.GridColor = [0.9 0.9 0.9];

% Disegniamo le barre raggruppate
hBar = bar(data, 'grouped', 'EdgeColor', 'none'); 

% Applichiamo i colori alle 2 colonne (Serie)
hBar(1).FaceColor = [0.20, 0.45, 0.65]; % Colonna 1: Delta V (Blu)
hBar(2).FaceColor = [0.85, 0.45, 0.25]; % Colonna 2: Delta T (Arancione)

% Configurazione asse X
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Strategy 1', 'Strategy 2', 'Strategy 3'}, ...
         'FontSize', 11); % Rimosso il $ del latex se usi testo normale, per una resa più pulita

% --- MODIFICA ETICHETTA ASSE Y ---
% Ora l'asse Y spiega entrambe le unità di misura convertite
ylabel('$\Delta V$ [km/s] \quad / \quad $\Delta T$ [$10^5$ s]', 'Interpreter', 'latex', 'FontSize', 12);
title('Comparison of Transfer Strategies', 'FontSize', 13, 'FontWeight', 'bold');

%% 8. LEGENDA (Sincronizzata con le colonne)
lgd = legend(hBar, {'$\Delta V$ [km/s]', '$\Delta T$ [days]'}, ...
    'Location', 'northeastoutside', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';

hold off;