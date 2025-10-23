function main_1_bSLIP_run
% Task 1 main – robust to either Signal Logging ('logsout') or To-Workspace (Timeseries).
clear; close all; clc;

setup_bSLIP_run;
mdl = 'bSLIP_run';
load_system(mdl);

% Make sure normal 5s run (Task 1)
flag_apex2apex = 0; %#ok<NASGU>
assignin('base','flag_apex2apex',flag_apex2apex);

% Solver (tight if needed)
set_param(mdl,'StopTime',num2str(t_end),'Solver','ode23t','MaxStep','1e-3','RelTol','1e-4');

simOut = sim(mdl,'ReturnWorkspaceOutputs','on');

% ---- fetch signals robustly ----
getV = make_fetch(simOut);

% time: use com_y time, else simOut.tout
t  = try_time(simOut, getV, {'com_y','y','dy','dcom_y'});
x  = try_data(getV, {'com_x','x'});
y  = try_data(getV, {'com_y','y'});

% feet (optional; fall back to NaNs if not logged)
fxL = try_data(getV, {'footL_x','foot_L_x','foot_x_L'}); 
fyL = try_data(getV, {'footL_y','foot_L_y','foot_y_L'});
fxR = try_data(getV, {'footR_x','foot_R_x','foot_x_R'});
fyR = try_data(getV, {'footR_y','foot_R_y','foot_y_R'});

% velocities (for energy if KE not logged)
dx = try_data(getV, {'dx','dcom_x','vx','com_vx'});
dy = try_data(getV, {'dy','dcom_y','vy','com_vy'});

% energies if available
KE     = try_data(getV, {'KE'});
PE_g   = try_data(getV, {'PE_g','PEgrav'});
PE_spr = try_data(getV, {'PE_spring','PEspr'});

% If energies not provided, compute basic ones
if any(isnan(KE)) || any(isnan(PE_g)) || any(isnan(PE_spr))
    % spring compression from leg lengths if available
    lL = try_data(getV, {'l_L','lL','legL_len'});
    lR = try_data(getV, {'l_R','lR','legR_len'});
    if any(isnan(lL)) || any(isnan(lR))
        compL = zeros(size(y));
        compR = zeros(size(y));
    else
        compL = max(0, l_nominal - lL);
        compR = max(0, l_nominal - lR);
    end

    if any(isnan(dx)) || any(isnan(dy))
        error('Could not find dx/dy to compute energies. Please log dx and dy (or dcom_x/dcom_y).');
    end

    KE     = 0.5*m*(dx.^2 + dy.^2);               % body only (feet small)
    PE_g   = m*g*y;                                % body only (good enough for HW1)
    PE_spr = 0.5*k_leg*(compL.^2 + compR.^2);      % sum of both springs
end
E_tot = KE + PE_g + PE_spr;

% ---------- Figure 1-1: spatial trajectories ----------
figure('Color','w'); hold on; box on;
plot(x, y, 'Color',[0 0 0],'LineWidth',2);
if ~all(isnan(fxL)), plot(fxL, fyL, 'Color',[.6 0 0],'LineWidth',2); end
if ~all(isnan(fxR)), plot(fxR, fyR, 'Color',[0 0 .6],'LineWidth',2); end
xlabel('x (m)'); ylabel('y (m)');
legend({'COM','foot L','foot R'},'Location','best');
title('Figure 1-1: Spatial trajectories');

% ---------- Figure 1-2: energies ----------
figure('Color','w'); hold on; box on;
plot(t, E_tot, 'Color',[0 0 0],'LineWidth',2);
plot(t, KE,    'Color',[0 0.45 0.74],'LineWidth',2);
plot(t, PE_g,  'Color',[0.85 0.33 0.10],'LineWidth',2);
plot(t, PE_spr,'Color',[0.47 0.67 0.19],'LineWidth',2);
xlabel('time (s)'); ylabel('Energy (J)');
legend('Total','Kinetic','Gravitational PE','Spring PE','Location','best');
title('Figure 1-2: Energy trajectories');

end

% ===== helpers =====
function getV = make_fetch(simOut)
    if any(strcmp(simOut.who,'logsout'))
        L = simOut.get('logsout');
        getV = @(name) fetch_from_logsout(L,name);
    else
        getV = @(name) fetch_from_ws(simOut,name);
    end
end
function ts = fetch_from_logsout(L,name)
    if any(strcmp({L.getElementNames},name))
        ts = L.getElement(name).Values;
    else
        ts = []; % not found
    end
end
function ts = fetch_from_ws(simOut,name)
    if any(strcmp(simOut.who,name))
        ts = simOut.get(name);
    else
        ts = [];
    end
end
function v = try_data(getV, names)
    v = nan;
    for i=1:numel(names)
        ts = getV(names{i});
        if ~isempty(ts)
            if isstruct(ts) && isfield(ts,'Time') && isfield(ts,'Data')
                v = ts.Data; return;
            elseif isa(ts,'timeseries')
                v = ts.Data; return;
            elseif isnumeric(ts)
                v = ts; return;
            end
        end
    end
    v = nan; v = v(ones(1,1)); % scalar NaN by default (will be promoted later)
end
function t = try_time(simOut, getV, preferNames)
    for i=1:numel(preferNames)
        ts = getV(preferNames{i});
        if ~isempty(ts) && isfield(ts,'Time')
            t = ts.Time; return;
        elseif isa(ts,'timeseries')
            t = ts.Time; return;
        end
    end
    if any(strcmp(simOut.who,'tout'))
        t = simOut.get('tout'); return;
    end
    error('Could not determine simulation time vector. Please log any signal (Timeseries) to get Time.');
end
