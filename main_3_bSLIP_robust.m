function main_3_bSLIP_robust

clear; close all; clc;

% model + setup
setup_bSLIP_run;
mdl = 'bSLIP_run'; load_system(mdl);
set_param(mdl,'ReturnWorkspaceOutputs','on');
%set_param(mdl,'EnableZeroCrossing','on');
set_param(mdl,'Solver','ode23t','MaxStep','1e-3','RelTol','1e-4');

% names to hunt
pref.dx = {'outdx','dx','dx_com','dcom_x'};
pref.dy = {'outdy','dy','dy_com','dcom_y'};
pref.y  = {'outcom_y','com_y','y','y_com'};

assignin('base','use_raibert',0); 
assignin('base','flag_apex2apex',1); 

% desired apex height + initial apex
y_target = 1.20;
phi_TD   = deg2rad(20);
y0=1.0; dy0=0; dx0=5.0;
assignin('base','y0',y0); assignin('base','dy0',dy0); assignin('base','dx0',dx0);

% estimate local slope once
dphi0 = deg2rad(1);
dH_dphi = estimate_dH_dphi(mdl,'setup_bSLIP_run',y0,dx0,phi_TD,dphi0,pref);

% step terrain between hops
base_h = evalin('base','plane_height');
steps  = [4 7 11 15];              
heights= [+0.05 -0.06 +0.04 -0.05];
h_total = base_h;

% hop loop
Tmax=10; hop_max=40; t_all=[]; y_all=[]; apex_log=[];

for hop=1:hop_max

    % terrain change?
    hit = find(steps==hop,1);
    if ~isempty(hit)
        h_total = h_total + heights(hit);
        assignin('base','plane_height',h_total);
    end

    assignin('base','phi_TD_cmd_rad',phi_TD);
    simOut = sim(mdl);

    [t,y,~] = get_ts_local(simOut,pref);
    t = t - t(1) + (t_all(end)+1e-6)*(~isempty(t_all));
    t_all = [t_all; t]; y_all = [y_all; y];
    y_apex = y(end); apex_log(end+1)=y_apex; 

    if t_all(end)>Tmax, break; end

    % deadbeat update for next hop
    phi_TD = deadbeat_phi_controller(phi_TD, y_apex, y_target, dH_dphi, [deg2rad(10) deg2rad(28)]);

    assignin('base','y0',y_apex);
    assignin('base','dy0',0);
    assignin('base','dx0', get_ts_local(simOut,pref,'dx','last'));
end

% plots
figure('Name','Task 3: apex vs hop'); clf; 
plot(0:numel(apex_log)-1, apex_log,'k.-'); hold on; yline(y_target,'r--');
xlabel('hop'); ylabel('apex height (m)'); box on;

figure('Name','Task 3: COM height (time)'); clf;
plot(t_all,y_all,'k-'); hold on; yline(y_target,'r--'); 
xlabel('time (s)'); ylabel('COM height (m)'); box on;
end

% helpers
function [t, y, dy] = get_ts_local(simOut, pref, which, pick)
if nargin<3, which=''; end, if nargin<4, pick=[]; end
try, yout = simOut.get('yout'); catch, yout = []; end
if ~isempty(yout)
    [y,t]  = grab_from_yout(yout, pref.y);
    [dy,~] = grab_from_yout(yout, pref.dy);
else
    try, logs = simOut.get('logsout'); catch, logs=[]; end
    if ~isempty(logs)
        [y,t]  = grab_from_logs(logs, pref.y);
        [dy,~] = grab_from_logs(logs, pref.dy);
    else
        [y,t]  = grab_from_ws(simOut, pref.y);
        [dy,~] = grab_from_ws(simOut, pref.dy);
    end
end
if isempty(which)
    return
else
    switch lower(which)
        case 'dx'
            v = try_one(simOut, pref.dx);
        case 'y'
            v = y;
        case 'dy'
            v = dy;
        otherwise
            error('unknown series request');
    end
    if isempty(pick), t = v; y = []; dy = []; return; end
    if ischar(pick) && strcmpi(pick,'last'), t=v(end); y=[]; dy=[]; else, t=v(pick); y=[]; dy=[]; end
end
end
function [v] = try_one(simOut, namelist)
v=[]; 
% yout
try, yout=simOut.get('yout'); for i=1:numel(yout.signals)
    lab=yout.signals(i).label; bn=yout.signals(i).blockName; dat=yout.signals(i).values.Data;
    if any(strcmpi(lab,namelist)) || any(endsWith(bn,namelist,'IgnoreCase',true)), v=dat; return; end
end, end
% logsout
try, logs=simOut.get('logsout');
    for i=1:numel(namelist), ts=logs.getElement(namelist{i}).Values; v=ts.Data; return; end
end
% ws
for i=1:numel(namelist)
    try, ts=simOut.get(namelist{i}); 
        if isa(ts,'timeseries'), v=ts.Data; return; end
        if isstruct(ts)&&isfield(ts,'Data'), v=ts.Data; return; end
    catch, end
end
end
function [v,t] = grab_from_yout(yout, namelist)
v=[]; t=[];
for i=1:numel(yout.signals)
    data=yout.signals(i).values; lab=yout.signals(i).label; bn=yout.signals(i).blockName;
    if any(strcmpi(lab,namelist)) || any(endsWith(bn,namelist,'IgnoreCase',true))
        t=data.Time; v=data.Data; return; end
end
end
function [v,t] = grab_from_logs(logs, namelist)
v=[]; t=[];
for i=1:numel(namelist)
    try, ts=logs.getElement(namelist{i}).Values; v=ts.Data; t=ts.Time; return; catch, end
end
end
function [v,t] = grab_from_ws(simOut, namelist)
v=[]; t=[];
for i=1:numel(namelist)
    try
        ts=simOut.get(namelist{i});
        if isa(ts,'timeseries'), v=ts.Data; t=ts.Time; return; end
        if isstruct(ts) && isfield(ts,'Time') && isfield(ts,'Data')
            v=ts.Data; t=ts.Time; return; 
        end
    catch, end
end
end
