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
hBar.CData(1, :) = [0.20, 0.45, 0.65]; % Colore per la Barra 1 (Blu)
hBar.CData(2, :) = [0.85, 0.45, 0.25]; % Colore per la Barra 2 (Arancione)

%% 7. CONFIGURAZIONE ASSI (2 Gruppi esatti)
set(gca, 'XTick', 1:2);

% Usiamo sprintf per inserire la variabile numerica 'th' all'interno della stringa LaTeX
etichetta_ottima = sprintf('$\\theta = %.2f^\\circ$', th); 

set(gca, 'XTickLabel', {'$\theta = 0^\circ$', etichetta_ottima}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

ylabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('Comparison of Transfer Strategies', 'FontSize', 13, 'FontWeight', 'bold');

%% 8. RISOLUZIONE DEL WARNING DELLA LEGENDA
% Creiamo due grafici invisibili (finti) sullo sfondo che hanno lo stesso 
% identico colore delle barre. Li useremo SOLO per ingannare la legenda.
dummy1 = bar(nan, 'FaceColor', [0.20, 0.45, 0.65], 'EdgeColor', 'none');
dummy2 = bar(nan, 'FaceColor', [0.85, 0.45, 0.25], 'EdgeColor', 'none');

% Ora passiamo alla legenda i due oggetti "dummy" specchio. Zero warning!
lgd = legend([dummy1, dummy2], {'Pericenter Burn', 'Optimal Burn'}, ...
    'Location', 'eastoutside');
lgd.FontSize = 10;
lgd.Box = 'on';

hold off;