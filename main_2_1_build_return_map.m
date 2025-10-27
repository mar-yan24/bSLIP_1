function main_2_1_build_return_map
clear; close all; clc;

% setup + model
setup_bSLIP_run;
mdl = 'bSLIP_run';
load_system(mdl);

% robust sim settings
set_param(mdl,'ReturnWorkspaceOutputs','on');
% set_param
set_param(mdl,'Solver','ode23t','MaxStep','1e-3','RelTol','1e-4');

% names we will look for 
pref.dx = {'outdx','dx','dx_com','dcom_x'};
pref.dy = {'outdy','dy','dy_com','dcom_y'};
pref.y  = {'outcom_y','com_y','y','y_com'};

% flags for model behaviors
assignin('base','use_raibert',0);
assignin('base','flag_apex2apex',1);

% reference energy: y_ref=1 m, dx_ref=5 m/s
y_ref = 1.0; dx_ref = 5.0;
E_ref = m*g*y_ref + 0.5*m*dx_ref^2;

% sweep settings
y_apex_vec      = (0.97:0.02:1.99).';   % column
phi_TD_list_deg = [15 20 25];

% --- Storage ---
RM = struct('phi_TD_deg',[],'y_i',[],'y_ip1',[]);
all_runs = [];

% loop over touchdown angles
for a = 1:numel(phi_TD_list_deg)

    phi_TD_cmd_rad = deg2rad(phi_TD_list_deg(a));
    assignin('base','phi_TD_cmd_rad',phi_TD_cmd_rad);

    y_i   = y_apex_vec;
    y_ip1 = nan(size(y_i));

    for k = 1:numel(y_apex_vec)

        % apex ICs for this run
        y0  = y_apex_vec(k);
        dy0 = 0;
        KE_x = max(E_ref - m*g*y0, 0);
        dx0  = sqrt(2*KE_x/m);

        assignin('base','y0', y0);
        assignin('base','dy0',dy0);
        assignin('base','dx0',dx0);

        % one hop; model should stop at the next apex
        simOut = sim(mdl);

        % robust fetch
        % --- Get signals robustly
        [t, y, dy] = get_ts_local(simOut, pref);   % use the local helper below
        
        % We *always* compute the next apex explicitly.
        % Skip a tiny window after t(1) (we start at an apex with dy=0; we want the *next* apex)
        j0 = find(t - t(1) > 1e-3, 1, 'first');           % 1 ms past start
        if isempty(j0), j0 = 2; end                       % safety
        
        s = sign(dy);
        idx_rel = find( s(j0-1:end-1) > 0 & s(j0:end) <= 0, 1, 'first');  % first +→– crossing
        assert(~isempty(idx_rel), 'No falling zero-crossing of dy found; apex not detected.');
        
        idx = (j0-1) + idx_rel;                           % absolute index of next apex
        
        y_i(k)   = y(1);
        y_ip1(k) = y(idx);


        % if model halted at apex, first sample = start apex, last = next apex
        if numel(t) >= 2 && t(end) > t(1)
            y_i(k)   = y(1);
            y_ip1(k) = y(end);
        else
            % fallback: detect first falling zero-crossing of dy
            idx = find_next_apex_idx_local(dy);
            assert(~isempty(idx),'No apex found (a=%d, k=%d).',a,k);
            y_i(k)   = y(1);
            y_ip1(k) = y(idx);
        end
    end

    RM(a).phi_TD_deg = phi_TD_list_deg(a);
    RM(a).y_i        = y_i;
    RM(a).y_ip1      = y_ip1;

    all_runs = [all_runs; table( repmat(phi_TD_list_deg(a),numel(y_i),1), y_i, y_ip1, ...
                                 'VariableNames',{'phi_TD_deg','y_i','y_ip1'})];
end

save('bSLIP_return_map.mat','RM','all_runs','y_apex_vec','phi_TD_list_deg','E_ref');
disp('Saved bSLIP_return_map.mat');
end

function [t, y, dy] = get_ts_local(simOut, pref)
% try outports, then logsout, then variables in simOut

try
    yout = simOut.get('yout');
catch, yout = []; end

if ~isempty(yout)
    [y,t]  = grab_from_yout(yout, pref.y);
    [dy,~] = grab_from_yout(yout, pref.dy);
    if isempty(t)
        error('Outports present but could not find y/dy. Check names.');
    end
    return
end

try
    logs = simOut.get('logsout');
catch, logs = []; end

if ~isempty(logs)
    [y,t]  = grab_from_logs(logs, pref.y);
    [dy,~] = grab_from_logs(logs, pref.dy);
    if ~isempty(t), return; end
end

[y,t]  = grab_from_ws(simOut, pref.y);
[dy,~] = grab_from_ws(simOut, pref.dy);
if isempty(t)
    error('Could not find y/dy in yout, logsout, or workspace variables.');
end
end

function [v,t] = grab_from_yout(yout, namelist)
v=[]; t=[];
for i=1:numel(yout.signals)
    data = yout.signals(i).values;
    lab  = yout.signals(i).label;
    bn   = yout.signals(i).blockName;
    if any(strcmpi(lab,namelist)) || any(endsWith(bn,namelist,'IgnoreCase',true))
        t = data.Time; v = data.Data; return;
    end
end
end

function [v,t] = grab_from_logs(logs, namelist)
v=[]; t=[];
for i=1:numel(namelist)
    try
        ts = logs.getElement(namelist{i}).Values;
        v = ts.Data; t = ts.Time; return;
    catch, end
end
end

function [v,t] = grab_from_ws(simOut, namelist)
v=[]; t=[];
for i=1:numel(namelist)
    try
        ts = simOut.get(namelist{i});
        if isa(ts,'timeseries'), v = ts.Data; t = ts.Time; return; end
        if isstruct(ts) && isfield(ts,'Time') && isfield(ts,'Data')
            v = ts.Data; t = ts.Time; return;
        end
    catch, end
end
end

function idx = find_next_apex_idx_local(dy)
s = sign(dy);
ii = find(s(1:end-1) > 0 & s(2:end) <= 0, 1, 'first');
if isempty(ii), idx = []; else, idx = ii+1; end
end
