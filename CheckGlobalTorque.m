function [residual_torque,residual_force_x,residual_force_y] = CheckGlobalTorque(xcm_stacked, ycm_stacked, Fx_stacked, Fy_stacked, Tz_stacked)
    % Inputs: 
    %   xcm_stacked, ycm_stacked : M x 1 vectors of particle CM positions
    %   Fx_stacked, Fy_stacked   : M x 1 vectors of baseline-subtracted net interaction forces
    %   Tz_stacked       : M x 1 vector of baseline-subtracted spin torques
    
    % 1. Orbital torques about global origin
    T_orbital = xcm_stacked .* Fy_stacked - ycm_stacked .* Fx_stacked;
    
    % 2. Global system torque sum
    residual_torque = sum(Tz_stacked) + sum(T_orbital);

    residual_force_x = sum(Fx_stacked);
    residual_force_y = sum(Fy_stacked);
    
    % 3. Output validation
    fprintf('Net System Interaction Force Residual: [%.2e, %.2e]\n', sum(Fx_stacked), sum(Fy_stacked));
    fprintf('Net System Global Yaw Torque Residual: %.2e\n', residual_torque);
end