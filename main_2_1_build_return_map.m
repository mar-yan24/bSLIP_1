clear; close all; clc;

setup_bSLIP_run;

% Apex detection parameters
apex_heights = 0.97:0.02:1.99;
phi_tgt_values = [15, 20, 25] * pi/180;

% Initialize storage
return_map_data = struct();

for p = 1:length(phi_tgt_values)
    phi_tgt = phi_tgt_values(p);
    y_apex_next = [];
    
    for i = 1:length(apex_heights)
        % Set initial apex state
        y0 = apex_heights(i);
        x0 = 0;
        dx0 = 5;
        dy0 = 0;
        
        % Configure for apex-to-apex
        flag_apex2apex = 1;
        
        % Run simulation
        try
            out = sim('bSLIP_run.slx', 'StopTime', '10');
            
            % Find next apex
            dy = out.dy_com;
            t = out.tout;
            y = out.y_com;
            
            % Detect apex crossings (dy changes from + to -)
            apex_idx = find(dy(1:end-1) > 0 & dy(2:end) <= 0);
            
            if length(apex_idx) >= 2
                y_apex_next(i) = y(apex_idx(2));
            else
                y_apex_next(i) = NaN;
            end
        catch
            y_apex_next(i) = NaN;
        end
    end
    
    return_map_data(p).phi_tgt = phi_tgt_values(p);
    return_map_data(p).y_apex = apex_heights;
    return_map_data(p).y_apex_next = y_apex_next;
end

save('bSLIP_return_map.mat', 'return_map_data');