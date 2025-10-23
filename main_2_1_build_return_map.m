function main_2_1_build_return_map
clear; close all; clc;

% === load Task-1 setup and model ===
setup_bSLIP_run;                 % your Task-1 setup (params & default ICs)
mdl = 'bSLIP_run';               % your model filename
load_system(mdl);

% --- Part 2 parameters from the rubric ---
y_apex_vec = 0.97:0.02:1.99;     % sweep (m)
phi_TD_list_deg = [15 20 25];    % touchdown hip angles
flag_apex2apex = 1;              % crucial: run apex-to-apex

% Reference energy (y_ref=1 m, dx_ref=5 m/s)
y_ref = 1.0; dx_ref = 5.0;
E_ref = m*g*y_ref + 0.5*m*dx_ref^2;

% Pre-allocate storage
RM = struct();   % Return-map data per phi_TD
all_runs = [];

for a = 1:numel(phi_TD_list_deg)
    phi_TD_deg = phi_TD_list_deg(a);          %#ok<NASGU>
    assignin('base','phi_TD_deg',phi_TD_list_deg(a));
    assignin('base','flag_apex2apex',flag_apex2apex);

    y_i  = y_apex_vec(:);
    y_ip1 = nan(size(y_i));
    apex0 = cell(size(y_i));
    apex1 = cell(size(y_i));

    for k = 1:numel(y_apex_vec)
        y0 = y_apex_vec(k);       %#ok<NASGU>
        dy0 = 0;                  %#ok<NASGU> % apex
        % pick dx0 to match reference energy
        KE_x = max(0, E_ref - m*g*y0);
        dx0 = sqrt( 2*KE_x/m );   %#ok<NASGU>
        % keep your existing horizontal x0; make sure both feet are initially off ground
        % (phi0_L,phi0_R,l0_L,l0_R, etc. may come from setup_bSLIP_run)
        % Push ICs to base:
        assignin('base','y0',y0);
        assignin('base','dy0',dy0);
        assignin('base','dx0',dx0);

        % Sim settings
        set_param(mdl,'StopTime','5','Solver','ode23t','MaxStep','1e-3','RelTol','1e-4');

        % Run
        simOut = sim(mdl,'ReturnWorkspaceOutputs','on');

        % Retrieve apex pairs:
        % Expect the model to log 'apex_y_pair' as [y_i, y_i1] the instant it stops
        % If you instead logged full time traces, compute next-apex here:
        if any(strcmp(simOut.who,'apex_y_pair'))
            ypair = simOut.get('apex_y_pair');
            y_i(k)   = ypair(1);
            y_ip1(k) = ypair(2);
        else
            % Fallback: infer next apex from logged y, dy
            L = getLogs(simOut);
            t  = L('y').Time;
            y  = L('y').Data; 
            dy = L('dy').Data;
            % find negative-going zero crossing of dy AFTER the start
            z = findZeroCrossNeg(dy,t);
            if numel(z)>=1
                y_i(k)   = y(1);
                y_ip1(k) = y(z(1));
            end
        end

        % Save states at apex (if you logged 'apex0' and 'apex1' structs)
        if any(strcmp(simOut.who','apex0'))
            apex0{k} = simOut.get('apex0');
        end
        if any(strcmp(simOut.who','apex1'))
            apex1{k} = simOut.get('apex1');
        end
    end

    RM(a).phi_TD_deg = phi_TD_list_deg(a);
    RM(a).y_i   = y_i;
    RM(a).y_ip1 = y_ip1;
    RM(a).apex0 = apex0;
    RM(a).apex1 = apex1;

    all_runs = [all_runs; table( repmat(phi_TD_list_deg(a),numel(y_i),1), y_i, y_ip1, ...
                                 'VariableNames',{'phi_TD_deg','y_i','y_ip1'})]; %#ok<AGROW>
end

save('bSLIP_return_map.mat','RM','all_runs','y_apex_vec','phi_TD_list_deg','E_ref');

disp('Saved bSLIP_return_map.mat');
end

% ---------- helpers ----------
function L = getLogs(simOut)
    % works with Signal Logging ('logsout') or To-Workspace timeseries
    if any(strcmp(simOut.who,'logsout'))
        ds = simOut.get('logsout');
        L = @(name) ds.getElement(name).Values;
    else
        L = @(name) simOut.get(name);
    end
end

function idx = findZeroCrossNeg(dy,t)
    s = sign(dy);
    dz = find(s(1:end-1)>0 & s(2:end)<=0);
    idx = dz+1;    % take the second sample around the crossing
end
