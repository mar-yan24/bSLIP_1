clear; close all; clc;
setup_bSLIP_run;  % loads params & ICs

mdl = 'bSLIP_run';
load_system(mdl);

% Set stop time
set_param(mdl,'StopTime',num2str(t_end));
% Tight solver if needed
set_param(mdl,'Solver','ode23t','MaxStep','1e-3','RelTol','1e-4');

% Push ICs to base (use your IC blocks or mask parameters in the model)
assignin('base','x0',x0);   assignin('base','y0',y0);
assignin('base','dx0',dx0); assignin('base','dy0',dy0);
assignin('base','phi0_L',phi0_L); assignin('base','phi0_R',phi0_R);
assignin('base','dphi0_L',dphi0_L); assignin('base','dphi0_R',dphi0_R);
assignin('base','l0_L',l0_L); assignin('base','l0_R',l0_R);
assignin('base','dl0_L',dl0_L); assignin('base','dl0_R',dl0_R);

% Run
simOut = sim(mdl,'ReturnWorkspaceOutputs','on');

% Grab logged data (replace names if yours differ)
t       = simOut.logsout.get('time').Values.Data;
com_x   = simOut.logsout.get('com_x').Values.Data;
com_y   = simOut.logsout.get('com_y').Values.Data;
footL_x = simOut.logsout.get('footL_x').Values.Data;
footL_y = simOut.logsout.get('footL_y').Values.Data;
footR_x = simOut.logsout.get('footR_x').Values.Data;
footR_y = simOut.logsout.get('footR_y').Values.Data;

KE      = simOut.logsout.get('KE').Values.Data;        % 0.5*m*(dx^2+dy^2) + foot inertias if you included them
PE_g    = simOut.logsout.get('PE_g').Values.Data;      % m*g*com_y (+ foot masses if modeled separately)
PE_spr  = simOut.logsout.get('PE_spring').Values.Data; % 0.5*k_leg*(max(0, l_nominal-l))^2 (sum for both legs)
E_tot   = KE + PE_g + PE_spr;

%% Figure 1-1: Spatial trajectories
figure('Color','w'); hold on; box on;
plot(com_x,   com_y,   'Color',[0 0 0],     'LineWidth',2);
plot(footL_x, footL_y, 'Color',[.6 0 0],    'LineWidth',2);
plot(footR_x, footR_y, 'Color',[0 0 .6],    'LineWidth',2);
xlabel('x (m)'); ylabel('y (m)');
legend('COM','foot L','foot R','Location','best'); title('bSLIP running trajectories');

%% Figure 1-2: Energies
figure('Color','w'); hold on; box on;
plot(t, E_tot, 'Color',[0 0 0],'LineWidth',2);
plot(t, KE,    'Color',[0 0.45 0.74],'LineWidth',2);
plot(t, PE_g,  'Color',[0.85 0.33 0.10],'LineWidth',2);
plot(t, PE_spr,'Color',[0.47 0.67 0.19],'LineWidth',2);
xlabel('time (s)'); ylabel('Energy (J)');
legend('Total','Kinetic','Gravitational PE','Spring PE','Location','best');
title('Energy trajectories');