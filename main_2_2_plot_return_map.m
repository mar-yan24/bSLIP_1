clear; close all; clc;

load('bSLIP_return_map.mat');

figure;
hold on;

% Unity line
plot([0.9, 2.0], [0.9, 2.0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);

% Return maps
styles = {'--', ':', '-.'};
for i = 1:length(return_map_data)
    valid = ~isnan(return_map_data(i).y_apex_next);
    plot(return_map_data(i).y_apex(valid), ...
         return_map_data(i).y_apex_next(valid), ...
         'LineStyle', styles{i}, 'Color', [0, 0, 0], 'LineWidth', 1);
end

xlabel('y_{apex,i} [m]');
ylabel('y_{apex,i+1} [m]');
legend('Unity', '\phi_{tgt} = 15°', '\phi_{tgt} = 20°', '\phi_{tgt} = 25°');
title('Apex Height Return Map');
grid on;
axis equal;
xlim([0.95, 2.0]);
ylim([0.95, 2.0]);