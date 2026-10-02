function [x_proj,y_proj,H,X_rot,Y_rot,Z_rot,N_surf_rot,dA_surf_rot,X_top,Y_top,Z_top,xcv_rot,ycv_rot,zcv_rot] = RotateAndShift(X_grid,Y_grid,Z_grid,x_nodes,y_nodes,N_surf,dA_surf,alpha,beta,gamma,f,xcm,ycm,zcm_height,to_plot_init,zcm_dist,xcog,ycog,thickness,is_flipped,xcv,ycv,zcv)

C = [xcm; ycm; zcm_height];
if ~is_flipped
    z_nodes = f(x_nodes,y_nodes,xcog,ycog) + zcm_height - zcm_dist; % z nodes of the BOTTOM contact line
else
    z_nodes = f(x_nodes,y_nodes,xcog,ycog) + zcm_height + zcm_dist; % z nodes of the TOP contact line
end

Z_grid = Z_grid + zcm_height - zcm_dist; % actual BOTTOM surface height
P = [x_nodes; y_nodes; z_nodes];

% Rx defines the rotation matrix around the X-axis by beta (roll)
Rx = [1, 0, 0; ...
      0, cos(beta), -sin(beta); ...
      0, sin(beta), cos(beta)];

% Ry defines the rotation matrix around the Y-axis by alpha (pitch)
Ry = [cos(alpha), 0, sin(alpha); ...
      0, 1, 0; ...
      -sin(alpha), 0, cos(alpha)];

% yaw
Rz = [cos(gamma), -sin(gamma), 0; ...
      sin(gamma),  cos(gamma), 0; ...
      0,           0,          1];

% Combined rotation matrix (yaw, then pitch, then roll)
R = Rx * Ry * Rz;

% project the body onto the plane
P_rotated = R * (P - C) + C;
% Extract the projected points
x_proj = P_rotated(1, :);
y_proj = P_rotated(2, :);
H = P_rotated(3,:)';

% now do bottom surface
% Store the original grid dimensions so we can reshape later
[rows, cols] = size(X_grid);
% point matrix
P = [X_grid(:)'; 
     Y_grid(:)'; 
     Z_grid(:)'];
% rotate
P_rotated = R * (P - C) + C;

% Extract rows and fold them back into 2D grids
X_rot = reshape(P_rotated(1, :), rows, cols);
Y_rot = reshape(P_rotated(2, :), rows, cols);
Z_rot = reshape(P_rotated(3, :), rows, cols);


% repeat for the top surface
P = [X_grid(:)'; 
     Y_grid(:)'; 
     Z_grid(:)' + thickness];
% rotate
P_rotated = R * (P - C) + C;

% Extract rows and fold them back into 2D grids
X_top = reshape(P_rotated(1, :), rows, cols);
Y_top = reshape(P_rotated(2, :), rows, cols);
Z_top = reshape(P_rotated(3, :), rows, cols);

% finally, repeat for center of volume
P = [xcv; ycv; zcv + zcm_height - zcm_dist];
P_rotated = R * (P - C) + C;
xcv_rot = P_rotated(1);
ycv_rot = P_rotated(2);
zcv_rot = P_rotated(3);


N_surf_rot = R * N_surf;
dA_surf_rot = R * dA_surf;




if to_plot_init
    figure;
    surf(X_rot, Y_rot, Z_rot, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on
    surf(X_top,Y_top,Z_top, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    scatter3(xcm,ycm,zcm_height,200,'filled','red');
    scatter3(xcv_rot,ycv_rot,zcv_rot,200,'filled','blue');
    scatter3(x_proj,y_proj,H);
    title('Rotated Object and Contact Line');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on; view(3);
end
end