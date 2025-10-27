function main_2_2_plot_return_map
% Load bSLIP_return_map.mat and draw the rubric's Figure 2-1.

clear; close all; clc;
S = load('bSLIP_return_map.mat');

figure('Name','Figure 2-1: Apex-height return map'); hold on; box on;
% Slope-1 line
ymin = min(S.y_apex_vec); ymax = max(S.y_apex_vec);
plot([ymin ymax],[ymin ymax],'Color',[.5 .5 .5],'LineWidth',1,'DisplayName','slope 1');

% Curves in required styles
styles = {'--',':','-.'};
for i = 1:numel(S.RM)
    plot(S.RM(i).y_i, S.RM(i).y_ip1, styles{i}, 'Color',[0 0 0], 'LineWidth',1, ...
        'DisplayName', sprintf('\\phi_{TD}=%d^\\circ', S.RM(i).phi_TD_deg));
end

xlabel('apex height y_i (m)'); ylabel('apex height y_{i+1} (m)');
legend('Location','southeast'); axis tight;
end
