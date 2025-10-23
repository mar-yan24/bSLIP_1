function main_2_2_plot_return_map
clear; close all; clc;

S = load('bSLIP_return_map.mat'); 
RM = S.RM;

figure('Color','w'); hold on; box on;

% slope-1 line covering the y-range in the file
ymin = min( cellfun(@(v)min(v), {RM.y_i}) );
ymax = max( cellfun(@(v)max(v), {RM.y_i}) );
plot([ymin ymax],[ymin ymax],'Color',[.5 .5 .5],'LineWidth',1);

% three maps: 15, 20, 25 deg
cols = [0 0 0; 0 .45 .74; .85 .33 .10];  % choose three distinct
markers = {'o','s','^'};

for k = 1:numel(RM)
    plot(RM(k).y_i, RM(k).y_ip1, markers{k}, ...
        'Color', cols(k,:), 'MarkerFaceColor','none','LineWidth',1.5);
end

xlabel('apex height y_i (m)');
ylabel('apex height y_{i+1} (m)');
legend({'slope 1','\phi_{TD}=15°','\phi_{TD}=20°','\phi_{TD}=25°'},'Location','best');
title('Figure 2-1: Apex-height return map');
grid on;

end
