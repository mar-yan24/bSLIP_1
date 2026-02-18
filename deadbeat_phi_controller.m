function phi_TD_new = deadbeat_phi_controller(phi_TD_current, y_apex, y_target, dH_dphi_est, phi_min_max)

    if nargin < 5, phi_min_max = [deg2rad(10) deg2rad(28)]; end

    err  = y_target - y_apex;
    if abs(dH_dphi_est) < 1e-6
        warning('Deadbeat: too small');
        dphi = 0;
    else
        dphi = err / dH_dphi_est;
    end

    phi_TD_new = phi_TD_current + dphi;
    phi_TD_new = min(max(phi_TD_new, phi_min_max(1)), phi_min_max(2));
end
