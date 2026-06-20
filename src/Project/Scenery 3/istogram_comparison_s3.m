scenery3_3D;
clc;
close all;

DV_pericenter = DV_pericentro;
DV_optimal = min_DV;
th = rad2deg(ottimo_th_i);

clearvars -except DV_pericenter DV_optimal th

% Incolonniamo i dati usando il punto e virgola ';'. 
% In questo modo creiamo 2 righe distinte (2 gruppi sull'asse X)
data = [
    DV_pericenter(1); 
    DV_optimal(1)
];

figure('Color', [1 1 1]); 
hold on; 
grid on;
ax = gca; 
ax.GridColor = [0.9 0.9 0.9];

% Disegniamo le barre. Avendo una sola colonna nella matrice, 
% hBar sarà un singolo oggetto che controlla TUTTE le barre del grafico.
hBar = bar(data, 'EdgeColor', 'none'); 

% --- COLORAZIONE INDIVIDUALE DELLE BARRE ---
% Poiché hBar è un oggetto unico, per colorare le barre in modo differenziato 
% dobbiamo attivare la proprietà 'CData' (Color Data) faccia per faccia.
hBar.FaceColor = 'flat'; 
hBar.CData(1, :) = [1.0, 0.1, 0.1] ; % Colore per l'iniezione al pericentro (Rosso)
hBar.CData(2, :) = [0.0, 0.4, 1.0]; % Colore per l'iniezione nel punto ottimale (Blu)

%% 7. CONFIGURAZIONE ASSI (2 Gruppi esatti)
set(gca, 'XTick', 1:2);

% Usiamo sprintf per inserire la variabile numerica 'th' all'interno della stringa LaTeX
etichetta_ottima = sprintf('$\\theta = %.2f^\\circ$', th); 

set(gca, 'XTickLabel', {'$\theta = 0^\circ $', etichetta_ottima}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

ylabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('$\Delta V$ Comparison for Injection Point in Scenario 3', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

%% 8. RISOLUZIONE DEL WARNING DELLA LEGENDA
% Creiamo due grafici invisibili (finti) sullo sfondo che hanno lo stesso 
% identico colore delle barre. Li useremo SOLO per ingannare la legenda.
dummy1 = bar(nan, 'FaceColor', [1.0, 0.1, 0.1], 'EdgeColor', 'none');
dummy2 = bar(nan, 'FaceColor', [0.0, 0.4, 1.0], 'EdgeColor', 'none');

% Ora passiamo alla legenda i due oggetti "dummy" specchio. Zero warning!
lgd = legend([dummy1, dummy2], {'Pericenter Burn', 'Optimal Burn'}, ...
    'Location', 'northeastoutside');
lgd.FontSize = 10;
lgd.Box = 'on';

hold off;