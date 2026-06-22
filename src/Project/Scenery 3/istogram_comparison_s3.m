scenery3_3D;
clc;
close all;

%estraggo solo i dati finali dallo scenario 3 nel caso 3D
DV_pericenter = DV_pericentro;
DV_optimal = min_DV;
th = rad2deg(ottimo_th_i);

clearvars -except DV_pericenter DV_optimal th

%Creo il dataset
data = [
    DV_pericenter(1); 
    DV_optimal(1)
];

%% GRAFICO 1: CONFRONTO PUNTO OTTIMO E PERICENTRO

figure('Color', [1 1 1]); 
hold on; 


hBar = bar(data, 'EdgeColor', 'none'); 



%modifico l'attributo CData per colorare le singole colonne di uno stesso raggruppamento
hBar.FaceColor = 'flat'; 
hBar.CData(1, :) = [1.0, 0.1, 0.1] ; % Colore per l'iniezione al pericentro (Rosso)
hBar.CData(2, :) = [0.0, 0.4, 1.0]; % Colore per l'iniezione nel punto ottimale (Blu)

set(gca, 'XTick', 1:2);

% Uso sprintf per creare la stringa contenente il theta ottimo
etichetta_ottima = sprintf('$\\theta = %.2f^\\circ$', th); 

set(gca, 'XTickLabel', {'$\theta = 0^\circ $', etichetta_ottima}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

ylabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('$\Delta V$ Comparison for Injection Point in Scenario 3', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

% Creiamo due colonne invisibili da associare alla legenda (altrimenti la legenda vede solo una 'bar' e viene mostrato solo il pericenter burn)
dummy1 = bar(nan, 'FaceColor', [1.0, 0.1, 0.1], 'EdgeColor', 'none');
dummy2 = bar(nan, 'FaceColor', [0.0, 0.4, 1.0], 'EdgeColor', 'none');

lgd = legend([dummy1, dummy2], {'Pericenter Burn', 'Optimal Burn'}, ...
    'Location', 'northeastoutside');
lgd.FontSize = 10;
lgd.Box = 'on';

hold off;