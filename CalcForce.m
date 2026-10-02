function [Fz,Fx,Fy,Tx_buoyancy,Ty_buoyancy,Tz_buoyancy,...
    Fx_buoyancy,Fy_buoyancy,Fz_buoyancy,...
    Fx_meniscus, Fy_meniscus,Fz_meniscus] = CalcForce(X_rot,Y_rot,Z_rot,valid_mask,...
    N_surf_rot,dA_surf_rot,mu,m,g,sigma,rho,w,xcm,ycm,zcm_height,nx,ny,is_flipped,thickness,...
    H_coll,r_x,r_y,lc,hs,Fx_baseline,Fy_baseline,gamma_rel,X_top,Y_top,Z_top,V_particle,xcv_rot,ycv_rot)

F_gravity = -m*g;

X_fin = X_rot(valid_mask);
Y_fin = Y_rot(valid_mask);
Z_fin = Z_rot(valid_mask);

X_top_fin = X_top(valid_mask);
Y_top_fin = Y_top(valid_mask);
Z_top_fin = Z_top(valid_mask);

nx_surf_rot = N_surf_rot(1,:)';
ny_surf_rot = N_surf_rot(2,:)';
nz_surf_rot = N_surf_rot(3,:)';

dAx_rot = dA_surf_rot(1,:)';
dAy_rot = dA_surf_rot(2,:)';
dAz_rot = dA_surf_rot(3,:)';




if ~is_flipped || thickness == 0      % contact line pinned at bottom OR infinitely thin


    Fz_buoyancy = -rho * g * sum(Z_fin .* dAz_rot);
    Fx_buoyancy = -rho * g * sum(Z_fin .* dAx_rot);
    Fy_buoyancy = -rho * g * sum(Z_fin .* dAy_rot);


    % Position relative to center of mass
    rx = X_fin - xcm;
    ry = Y_fin - ycm;
    rz = Z_fin - zcm_height;

    % Buoyancy torque: r x dF
    Tx_buoyancy = -rho * g * sum(Z_fin .* ...
        (ry .* dAz_rot - rz .* dAy_rot));

    Ty_buoyancy = -rho * g * sum(Z_fin .* ...
        (rz .* dAx_rot - rx .* dAz_rot));

    Tz_buoyancy = -rho * g * sum(Z_fin .* ...
        (rx .* dAy_rot - ry .* dAx_rot));


else                % contact line pinned at top
    % Entire bottom + side surface is wetted.
    % Top surface is NOT wetted.
    %
    % dA_rot is the vector-area element for the bottom surface.
    % The top surface has the same vector-area elements because it is a
    % node-for-node translation of the bottom surface, but they flip so as
    % to still point inward
    
    dAx_top = -dAx_rot;
    dAy_top = -dAy_rot;
    dAz_top = -dAz_rot;
    
    
    % Pressure force on the UNWETTED top surface
    
    Fx_top = -rho * g * sum(Z_top_fin .* dAx_top);
    Fy_top = -rho * g * sum(Z_top_fin .* dAy_top);
    Fz_top = -rho * g * sum(Z_top_fin .* dAz_top);
    
    
    % Buoyancy force on the wetted surface
    %
    % complete closed particle:
    %     F_closed = [0, 0, rho*g*V_particle]
    %
    % wetted surface = closed surface - top surface

    Fx_closed = 0;
    Fy_closed = 0;
    Fz_closed = rho * g * V_particle;
    
    Fx_buoyancy = Fx_closed - Fx_top;
    Fy_buoyancy = Fy_closed - Fy_top;
    Fz_buoyancy = Fz_closed - Fz_top;
    
    
    % Position relative to center of mass
    
    rx_top = X_top_fin - xcm;
    ry_top = Y_top_fin - ycm;
    rz_top = Z_top_fin - zcm_height;
    
    
    % Hydrostatic torque on the top surface
    
    Tx_top = -rho * g * sum(Z_top_fin .* ...
        (ry_top .* dAz_top - rz_top .* dAy_top));
    
    Ty_top = -rho * g * sum(Z_top_fin .* ...
        (rz_top .* dAx_top - rx_top .* dAz_top));
    
    Tz_top = -rho * g * sum(Z_top_fin .* ...
        (rx_top .* dAy_top - ry_top .* dAx_top));
    
    
    % Buoyancy torque

    % Lever arm from center of mass to center of volume
    rx_cv = xcv_rot - xcm;
    ry_cv = ycv_rot - ycm;

    Tx_closed = ry_cv * Fz_closed;
    Ty_closed = -rx_cv * Fz_closed;
    Tz_closed = 0;

    Tx_buoyancy = Tx_closed - Tx_top;
    Ty_buoyancy = Ty_closed - Ty_top;
    Tz_buoyancy = Tz_closed - Tz_top;


end


% normals oriented opposite of in the BEM solve
nfx = -nx;
nfy = -ny;


% Force all arrays to be column vectors to fix dimension mismatches
nfx = nfx(:);
nfy = nfy(:);
mu = mu(:);
hs = hs(:);
w = w(:);

% tangent vectors
tfx = -nfy;
tfy = nfx;

% Normal slope (from BEM solution)
hn = -mu;

% Tangential slope hs = dh/ds - passed in

% O(h^2) capillary forces (linearized BEM) -> use the quadratic expansion

Fz_meniscus = sigma * sum( hn .* w );
Fx_meniscus = sigma * sum( ( (1 + 0.5*(hs.^2 - hn.^2)) .* nfx - (hn .* hs) .* tfx ) .* w );
Fy_meniscus = sigma * sum( ( (1 + 0.5*(hs.^2 - hn.^2)) .* nfy - (hn .* hs) .* tfy ) .* w );


% baseline force subtraction

% rotate by relative angle of this particle to particle 1
Rz = [cos(gamma_rel), -sin(gamma_rel); ...
      sin(gamma_rel),  cos(gamma_rel)];

F_sub = Rz * [Fx_baseline; Fy_baseline];

Fx_meniscus = Fx_meniscus - F_sub(1);
Fy_meniscus = Fy_meniscus - F_sub(2);



Fz = F_gravity + Fz_buoyancy + Fz_meniscus;
Fx = Fx_buoyancy + Fx_meniscus;
Fy = Fy_buoyancy + Fy_meniscus;



end

