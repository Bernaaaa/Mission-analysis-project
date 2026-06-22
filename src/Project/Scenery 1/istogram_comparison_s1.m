clear; clc; close all;
%% DATA FROM STRATEGY 1: STANDARD
scenery1_standard; 
clc
close all

DV_pl_standard = abs( DeltaVP_BLT);
DV_BT_standard = abs(DeltaV1BT) + abs(DeltaV2BT); 
DV_cp_standard = abs(DeltaVCP);
DV_total_standard = abs(DeltaV_total);
DT_standard = DeltaT_TOT;

% creating temp files to store variables
save('temp_strat1.mat', 'DV_pl_standard', 'DV_BT_standard', 'DV_cp_standard', 'DV_total_standard', 'DT_standard');
clear
%% DATA FROM STRATEGY 2: DELTA T
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
%% DATA FROM STRATEGY 3: DELTA V
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
%% COSTRUZIONE DATA SET
load('temp_strat1.mat');
load('temp_strat2.mat');
load('temp_strat3.mat')

% Deleting temp files
delete('temp_strat1.mat');
delete('temp_strat2.mat');
delete('temp_strat3.mat');

colonna_Standard = [DV_total_standard(1); DV_cp_standard(1); DV_pl_standard(1); DV_BT_standard(1)];
colonna_DeltaT = [DV_total_Delta_T(1); DV_cp_Delta_T(1); DV_pl_Delta_T(1); DV_BT_Delta_T(1)];
colonna_DeltaV = [DV_total_deltaV(1); DV_cp_deltaV(1); DV_pl_deltaV(1);DV_BT_deltaV(1)];

data = [colonna_Standard, colonna_DeltaT, colonna_DeltaV];
%% ISTOGRAMMA 1: DV MANEUVER COMPARISON

figure('Color', [1 1 1]); 
hold on; 
grid on;

blu      = [31, 119, 180] / 255;
arancio  = [255, 127, 14] / 255;
verde    = [44, 160, 44]  / 255;

hBar = barh(data, 'grouped', 'EdgeColor', 'none'); 

%colour palette
hBar(1).FaceColor = arancio; % Strategy 1 (Orange)
hBar(2).FaceColor = blu;     % Strategy 2 (Blu)
hBar(3).FaceColor = verde;   % Srategy 3 (Green)

%Creating Y-axis
set(gca, 'YTick', 1:4);
set(gca, 'YTickLabel', {'$\Delta V_{TOT}$', '$\Delta V_{per}$', '$\Delta V_{plane}$', '$\Delta V_{BT}$'}, ...
         'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Setting up labels and immage title
xlabel('$\Delta V$ [km/s]', 'Interpreter', 'latex', 'FontSize', 12);
title('$\Delta V$ Strategy Comparison', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');

%Creating legend
lgd = legend(hBar, {'First Strategy', 'Second Strategy', 'Third Strategy'}, ...
    'Location', 'northeastoutside', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.String = 'Transfer Strategies';
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';

hold off;
%% ISTOGRAMMA 2: STRATEGY DV E DT COMPARISON
%rearranging data
DV_data = [DV_total_standard(1); DV_total_Delta_T(1); DV_total_deltaV(1)];
DT_data = [DT_standard(1); DT_Delta_T(1); DT_DeltaV(1)];

data = [DV_data, DT_data];
%scaling time values
data(:, 2) = data(:, 2) ./ 10^5; 


figure('Color', [1 1 1]); 
hold on; 

%colour palette:
hBar = bar(data, 'grouped', 'EdgeColor', 'none'); 
hBar(1).FaceColor = [0.0, 0.4, 1.0]; % Column 1: Delta V (Blu)
hBar(2).FaceColor = [1.0, 0.1, 0.1]; % Column 2 : Delta T (Red)

% Creating X-axis
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Strategy 1', 'Strategy 2', 'Strategy 3'}, ...
         'FontSize', 11); 

% Creating y-label
ylabel('$\Delta V$ [km/s] \quad / \quad $\Delta T$ [$10^5$s]', 'Interpreter', 'latex', 'FontSize', 12);

%Setting up legend
lgd = legend(hBar, {'$\Delta V$', '$\Delta T$'}, ...
    'Location', 'northwest', 'Interpreter', 'latex');
lgd.FontSize = 10;
lgd.Title.Interpreter = 'latex';
lgd.Box = 'on';

hold off;