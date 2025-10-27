function dH = estimate_dH_dphi(mdl, setupFcn, y_apex, dx_apex, phi0, dphi, pref)
% Estimate dH/dphi around phi0 from two apex->apex runs.

    feval(setupFcn);
    load_system(mdl);
    set_param(mdl,'ReturnWorkspaceOutputs','on');
    %set_param(mdl,'EnableZeroCrossing','on');
    assignin('base','use_raibert',0);
    assignin('base','flag_apex2apex',1);

    % Common apex ICs
    assignin('base','y0',y_apex);
    assignin('base','dy0',0);
    assignin('base','dx0',dx_apex);

    % minus
    assignin('base','phi_TD_cmd_rad',phi0 - dphi);
    simOutM = sim(mdl);
    [~, yM] = get_ts_local(simOutM, pref);   % local helper below
    y_next_m = yM(end);

    % plus
    assignin('base','phi_TD_cmd_rad',phi0 + dphi);
    simOutP = sim(mdl);
    [~, yP] = get_ts_local(simOutP, pref);
    y_next_p = yP(end);

    dH = (y_next_p - y_next_m) / (2*dphi);
end

% local helper (identical to the one in main_2_1 for convenience)
function [t, y, dy] = get_ts_local(simOut, pref)
try, yout = simOut.get('yout'); catch, yout = []; end
if ~isempty(yout)
    [y,t]  = grab_from_yout(yout, pref.y);
    [dy,~] = grab_from_yout(yout, pref.dy); return;
end
try, logs = simOut.get('logsout'); catch, logs=[]; end
if ~isempty(logs)
    [y,t]  = grab_from_logs(logs, pref.y);
    [dy,~] = grab_from_logs(logs, pref.dy); if ~isempty(t), return; end
end
[y,t]  = grab_from_ws(simOut, pref.y);
[dy,~] = grab_from_ws(simOut, pref.dy);
if isempty(t), error('get_ts_local: could not find y/dy'); end
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
        if isstruct(ts) && isfield(ts,'Time'), v=ts.Data; t=ts.Time; return; end
    catch, end
end
end
