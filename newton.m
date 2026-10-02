clear
clc;
close all

%% Initialize parameters. All units cm g s

% fluid parameters and constants
sigma = 72;               % surface tension
rho = 1;                  % density
g = 981;                  % gravity
lc = sqrt(sigma/(rho*g)); % capillary length
gamma_E = 0.5772156649;   % Euler-Mascheroni constant

% for iteration: tolerance for stopping, max iterations
tol = 1e-9;
max_iters = 20; % should not need nearly this many, but just to be safe

% step sizes for centered differences in Jacobian
dz = tol;
dalpha = tol;
dbeta = tol;

limits = [-0.0075 0.0075]; % colorbar limits for plotting fluid height
plotCM = 1; % toggle to plot the CoM of each particle (1 = plot, 0 = don't)

Nmin = 150; % minimum number of boundary points per particle
Nmax = 400; % maximum number of boundary points per particle
C_domain = 0.01; % domain grid spacing is 1% of the characteristic size of the particle assembly    
C_BEM = 0.05; % want 1/C_BEM = 20 boundary points across the particle/gap length scale
C_curv = 0.05; % want 1/C_curv = 20 boundary points per min radius of curvature

view_init = 0; % flag to plot initial meniscus solve before equilibrium is found
to_plot_init = 0; % flag to plot initial rotation stuff. Keep this turned off if doing full iteration! This is for testing
check_manufactured = 1; % check BEM solve against a manufactured solution

%% GUI stuff

% set particle number, mass and thickness, choose pinning
[m, thickness, M, is_flipped, snapToGrid, gridSpacing, inputMethod] = getParametersGUI();

% set particle positions (centers of geometry) and characteristic length
if strcmp(inputMethod,'graphical')
    [xcog_stacked, ycog_stacked, L] = particlePositionGraphicalGUI(gridSpacing, snapToGrid, M);
else
    [xcog_stacked, ycog_stacked, L] = particlePositionWriteInGUI(M);
end

% particle surface equation
dev = surfaceWavinessGUI(); % dev is the dimensionless amplitude of deformation relative to L: delta z / L = dev, where delta z is max amplitude
f = @(x,y,xcog,ycog) dev * L * (((x-xcog)/L).^2 - ((y-ycog)/L).^2);

% get contour and center of mass displacement
[x_nodes, y_nodes, xcm_dist, ycm_dist, zcm_dist, user_set_cm, N, perimeter, Rmin, ds] = drawContourGUI(L, snapToGrid, gridSpacing, M, xcog_stacked, ycog_stacked, dev, Nmin, Nmax, C_BEM, C_curv, thickness); % xcm_dist and ycm_dist are distance of the CoM to the CoG (how displaced is the CoM?)

% zcm_dist is the vertical distance FROM the surface, z = f(xcm,ycm), TO the center of mass
if user_set_cm
    zcm_dist_stacked = zcm_dist*ones(M,1);
else
    zcm_dist_stacked = zeros(M,1); % if not, we'll compute the CoV later
end

% bookkeeping of stacked representations of particles and their CoMs
m_stacked = m*ones(M,1);
x_nodes = x_nodes';
y_nodes = y_nodes';
x_nodes_stacked = repmat(x_nodes,M,1) + xcog_stacked;
y_nodes_stacked = repmat(y_nodes,M,1) + ycog_stacked;

xcm_dist_stacked = xcm_dist*ones(M,1);
ycm_dist_stacked = ycm_dist*ones(M,1);
xcm_stacked = xcog_stacked + xcm_dist_stacked;
ycm_stacked = ycog_stacked + ycm_dist_stacked;

% for domain plotting
[min_x,max_x] = bounds(x_nodes_stacked,'all');
[min_y,max_y] = bounds(y_nodes_stacked,'all');
char_scale = max(max_x-min_x,max_y-min_y); % for plotting
grid_range = 1.2*char_scale; 

% domain points: resolve overall particle-assembly scale
dx_domain = C_domain*char_scale;
N_domain = ceil(grid_range/dx_domain) + 1;

% user chooses heights and initial rotations
[zcm_height_stacked, alpha_stacked, beta_stacked, gamma_stacked] = initialConditionsGUI(M, x_nodes_stacked, y_nodes_stacked, xcog_stacked, ycog_stacked, xcm_stacked, ycm_stacked);

% convert to radians
alpha_stacked = alpha_stacked*pi/180;
beta_stacked = beta_stacked*pi/180;
gamma_stacked = gamma_stacked*pi/180;

gamma_rel_stacked = gamma_stacked - gamma_stacked(1); % relative yaw angle of all particles compared to particle 1


%% initialize

is_inside = cell(1,M);
is_inside_uniform = cell(1,M);
N_surf_rot_stacked = cell(1,M);
dA_surf_rot_stacked = cell(1,M);
N_surf_stacked = cell(1,M);
dA_surf_stacked = cell(1,M);
valid_mask_stacked = cell(1,M);

% centers of volume and volumes -- to be computed
xcv_stacked = zeros(M,1);
ycv_stacked = zeros(M,1);
zcv_stacked = zeros(M,1);
xcv_rot_stacked = zeros(M,1);
ycv_rot_stacked = zeros(M,1);
zcv_rot_stacked = zeros(M,1);
V_stacked = zeros(M,1);
    
x_proj_stacked = zeros(M,N);
y_proj_stacked = zeros(M,N);
H_stacked = zeros(N,M); % not a typo, has to be N then M for matrix multiplication later
X_grid_stacked = cell(1,M);
Y_grid_stacked = cell(1,M);
Z_grid_stacked = cell(1,M);

% these are grids that represent the BOTTOM surface
X_rot_stacked = cell(1,M);
Y_rot_stacked = cell(1,M);
Z_rot_stacked = cell(1,M);

% these represent the TOP surface
X_top_stacked = cell(1,M);
Y_top_stacked = cell(1,M);
Z_top_stacked = cell(1,M);


% Interpolate uniform area grid, get normals and area elements

for p=1:M % for each particle

    x_nodes = x_nodes_stacked(p,:);
    y_nodes = y_nodes_stacked(p,:);
    alpha = alpha_stacked(p);
    beta = beta_stacked(p);
    gamma = gamma_stacked(p);
    zcm_height = zcm_height_stacked(p);
    xcm = xcm_stacked(p);
    ycm = ycm_stacked(p);
    xcog = xcog_stacked(p);
    ycog = ycog_stacked(p);

    [X,Y] = meshgrid(x_nodes,y_nodes);
    is_inside{p} = inpolygon(X, Y, x_nodes, y_nodes);
    Z = f(X,Y,xcog,ycog);

    % apply mask before rotation
    X(~is_inside{p}) = NaN;
    Y(~is_inside{p}) = NaN;
    Z(~is_inside{p}) = NaN;

    % isolate valid unique points and build the interpolant
    valid_idx = ~isnan(Z);
    X_val = X(valid_idx);
    Y_val = Y(valid_idx);
    Z_val = Z(valid_idx);
    coords = [X_val(:), Y_val(:)];
    [~, unique_idx, ~] = unique(coords, 'rows');
    X_clean = X_val(unique_idx);
    Y_clean = Y_val(unique_idx);
    Z_clean = Z_val(unique_idx);

    F = scatteredInterpolant(X_clean, Y_clean, Z_clean, 'linear', 'none');

    % uniform flat grid
    grid_resolution = L / N;

    % Find the furthest distance from the CoM to any point on the object
    max_dist_x = max(abs(X_val - xcm));
    max_dist_y = max(abs(Y_val - ycm));

    % Calculate how many full steps are needed to cover that distance
    n_steps_x = ceil(max_dist_x / grid_resolution);
    n_steps_y = ceil(max_dist_y / grid_resolution);

    % Build a grid radiating outward from the CoM
    x_bounds = xcm + (-n_steps_x : n_steps_x) * grid_resolution;
    y_bounds = ycm + (-n_steps_y : n_steps_y) * grid_resolution;

    % uniform non-rotated grid; this will remain fixed
    [X_grid, Y_grid] = meshgrid(x_bounds, y_bounds);

    % evaluate the interpolant ON the uniform grid
    Z_grid = F(X_grid, Y_grid);
    is_inside_uniform{p} = inpolygon(X_grid, Y_grid, x_nodes, y_nodes);
    Z_grid(~is_inside_uniform{p}) = NaN;

    % compute surface slopes on the full grid
    [fx,fy] = gradient(Z_grid, grid_resolution, grid_resolution);
    valid_mask_stacked{p} = ~isnan(Z_grid) & ~isnan(fx) & ~isnan(fy);

    % vectors
    X_fin = X_grid(valid_mask_stacked{p});
    Y_fin = Y_grid(valid_mask_stacked{p});
    Z_fin = Z_grid(valid_mask_stacked{p});

    fx_fin = fx(valid_mask_stacked{p});
    fy_fin = fy(valid_mask_stacked{p});

    cell_area = grid_resolution^2;

    % Unit normal to surface
    nz_surf = 1 ./ sqrt(1 + fx_fin.^2 + fy_fin.^2);
    nx_surf = -fx_fin .* nz_surf;
    ny_surf = -fy_fin .* nz_surf;

    N_surf = [nx_surf(:)'; 
          ny_surf(:)'; 
          nz_surf(:)'];

    N_surf_stacked{p} = N_surf;

    % Actual vector surface area element
    dAx = -fx_fin * cell_area;
    dAy = -fy_fin * cell_area;
    dAz =  ones(size(fx_fin))*cell_area;

    dA_surf = [dAx(:)';
               dAy(:)';
               dAz(:)'];

    dA_surf_stacked{p} = dA_surf;

    X_grid_stacked{p} = X_grid;
    Y_grid_stacked{p} = Y_grid;
    Z_grid_stacked{p} = Z_grid;

    % Calculate volume and center of volume
    inside = valid_mask_stacked{p};
    dA = grid_resolution^2;
    % Footprint area
    A = nnz(inside) * dA;
    % Volume
    V_stacked(p) = thickness * A;
    % Center of volume - relative to the origin
    xcv_stacked(p) = sum(X_grid(inside)) * dA / A;
    ycv_stacked(p) = sum(Y_grid(inside)) * dA / A;
    zcv_stacked(p) = sum(Z_grid(inside)) * dA / A + thickness/2;
    
    % zcm_dist is the position of the ZCM relative to the BOTTOM SURFACE
    if ~user_set_cm % if the user didn't assign the CoM, CoM = CoV = CoG, so use the computed z-value
        zcm_dist_stacked(p) = zcv_stacked(p) - f(xcv_stacked(p),ycv_stacked(p),xcog,ycog);
    end
    zcm_dist = zcm_dist_stacked(p);

    % put bodies in current orientation
    [x_proj,y_proj,H,X_rot,Y_rot,Z_rot,N_surf_rot,dA_surf_rot,X_top,Y_top,Z_top,xcv_rot,ycv_rot,zcv_rot] = RotateAndShift(X_grid,Y_grid,...
                                                                                Z_grid,x_nodes,y_nodes,N_surf,...
                                                                                dA_surf,alpha,beta,gamma,f,xcm,ycm,...
                                                                                zcm_height,to_plot_init,zcm_dist,xcog,ycog,thickness,is_flipped,xcv_stacked(p),ycv_stacked(p),zcv_stacked(p));

    x_proj_stacked(p,:) = x_proj;
    y_proj_stacked(p,:) = y_proj;
    H_stacked(:,p) = H; % z-coordinates of the pinned contact line
    X_rot_stacked{p} = X_rot;
    Y_rot_stacked{p} = Y_rot;
    Z_rot_stacked{p} = Z_rot;
    N_surf_rot_stacked{p} = N_surf_rot;
    dA_surf_rot_stacked{p} = dA_surf_rot;
    X_top_stacked{p} = X_top;
    Y_top_stacked{p} = Y_top;
    Z_top_stacked{p} = Z_top;
    xcv_rot_stacked(p) = xcv_rot;
    ycv_rot_stacked(p) = ycv_rot;
    zcv_rot_stacked(p) = zcv_rot;


end


if thickness ~= 0
    fprintf('Particle density: %.2f g/cm^3\n',m_stacked(1)/V_stacked(1));
else
    fprintf('Particle density: N/A (0 thickness)\n');
end


fprintf('Forces [g cm / s^2] = dynes, torques [g cm^2 / s^2] = dyne-centimeters\n\n');



% error handling: contours cannot intersect or overlap. Pairwise check
overlap_matrix = false(M,M);
polys = cell(M,1);
for p = 1:M
    polys{p} = polyshape(x_proj_stacked(p,:), y_proj_stacked(p,:));
end

for p = 1:M
    for q = p+1:M
        overlap_matrix(p,q) = overlaps(polys{p},polys{q});
        overlap_matrix(q,p) = overlap_matrix(p,q);
    end
end
if any(overlap_matrix,'all')
    figure;
    scatter(x_proj_stacked,y_proj_stacked);
    hold on
    scatter(xcm_stacked,ycm_stacked);
    error('Particles overlap!');
end





%% Baseline spurious force for a single particle (k=1)

% make cell representations for just this single particle (clunky but oh well)
N_surf_rot_lone_cell = cell(1,1);
dA_surf_rot_lone_cell = cell(1,1);
N_surf_rot_lone_cell{1} = N_surf_rot_stacked{1};
dA_surf_rot_lone_cell{1} = dA_surf_rot_stacked{1};

N_surf_lone_cell = cell(1,1);
dA_surf_lone_cell = cell(1,1);
N_surf_lone_cell{1} = N_surf_stacked{1};
dA_surf_lone_cell{1} = dA_surf_stacked{1};
    
X_grid_lone_cell = cell(1,1);
Y_grid_lone_cell = cell(1,1);
Z_grid_lone_cell = cell(1,1);
X_grid_lone_cell{1} = X_grid_stacked{1};
Y_grid_lone_cell{1} = Y_grid_stacked{1};
Z_grid_lone_cell{1} = Z_grid_stacked{1};

% these are grids
X_rot_lone_cell = cell(1,1);
Y_rot_lone_cell = cell(1,1);
Z_rot_lone_cell = cell(1,1);
X_rot_lone_cell{1} = X_rot_stacked{1};
Y_rot_lone_cell{1} = Y_rot_stacked{1};
Z_rot_lone_cell{1} = Z_rot_stacked{1};

X_top_lone_cell = cell(1,1);
Y_top_lone_cell = cell(1,1);
Z_top_lone_cell = cell(1,1);
X_top_lone_cell{1} = X_top_stacked{1};
Y_top_lone_cell{1} = Y_top_stacked{1};
Z_top_lone_cell{1} = Z_top_stacked{1};

valid_mask_lone_cell = cell(1,1);
valid_mask_lone_cell{1} = valid_mask_stacked{1};

% BEM solve for one particle
[mu_lone,rx_lone,ry_lone,nx_lone,ny_lone,w_lone,H_coll_lone,hs_lone] = mu_solve(x_proj_stacked(1,:),y_proj_stacked(1,:),N,H_stacked(:,1),lc,gamma_E,1);


% calculate total force
[Fz,~,~,Tx_buoyancy,Ty_buoyancy,~,...
~,~,~,...
~, ~,~] = CalcForce(X_rot_lone_cell{1},Y_rot_lone_cell{1},Z_rot_lone_cell{1},valid_mask_lone_cell{1},...
N_surf_rot_lone_cell{1},dA_surf_rot_lone_cell{1},mu_lone,m,g,sigma,rho,w_lone,xcm_stacked(1),ycm_stacked(1),zcm_height_stacked(1),...
nx_lone,ny_lone,is_flipped,thickness,H_coll_lone,rx_lone,ry_lone,lc,hs_lone,0,0,0,...
X_top_lone_cell{1},Y_top_lone_cell{1},Z_top_lone_cell{1},V_stacked(1),xcv_rot_stacked(1),ycv_rot_stacked(1));

% calculate torques
[Tx_men,Ty_men,~] = CalcTorque(mu_lone,w_lone,sigma,rx_lone,ry_lone,xcm_stacked(1),ycm_stacked(1),nx_lone,ny_lone,...
    H_coll_lone,zcm_height_stacked(1),hs_lone,0);
Tx = Tx_men + Tx_buoyancy;
Ty = Ty_men + Ty_buoyancy;

Fz_baseline = Fz;
Tx_baseline = Tx;
Ty_baseline = Ty;



hWaitbar = waitbar(0,'Finding lone particle equilibrium...', 'Name', 'Status');


% compute Jacobian of lone particle prior to equilibrium
J_lone = Jacobian(1,x_nodes_stacked(1,:),y_nodes_stacked(1,:),xcm_stacked(1),ycm_stacked(1),...
    xcog_stacked(1,:),ycog_stacked(1,:),zcm_dist_stacked(1),gamma_stacked(1),m,zcm_height_stacked(1),...
    N_surf_lone_cell,dA_surf_lone_cell,X_grid_lone_cell,Y_grid_lone_cell,Z_grid_lone_cell,...
    alpha_stacked(1),beta_stacked(1),x_proj_stacked(1,:),y_proj_stacked(1,:),H_stacked(:,1),...
    X_rot_lone_cell,Y_rot_lone_cell,Z_rot_lone_cell,N_surf_rot_lone_cell,...
    dA_surf_rot_lone_cell,dz,dbeta,dalpha,f,to_plot_init,is_flipped,thickness,...
    N,lc,gamma_E,valid_mask_lone_cell,g,sigma,rho,0,0,0,0,xcv_stacked(1),ycv_stacked(1),zcv_stacked(1),...
    X_top_lone_cell,Y_top_lone_cell,Z_top_lone_cell,xcv_rot_stacked(1),ycv_rot_stacked(1),zcv_rot_stacked(1),V_stacked);


iters = 0;

zcm_height_lone = zcm_height_stacked(1);
beta_lone = beta_stacked(1);
alpha_lone = alpha_stacked(1);

while 1

    % Newton solve
    q = reshape([zcm_height_lone, beta_lone, alpha_lone].', [], 1);
    R = reshape([Fz_baseline, Tx_baseline, Ty_baseline].', [], 1);

    dq = -J_lone \ R;

    q = q + dq;

    zcm_height_lone = q(1:3:end);
    beta_lone = q(2:3:end);
    alpha_lone = q(3:3:end);


    % shift and rotate particle accorindgly

    xcm = xcm_stacked(1);
    ycm = ycm_stacked(1);
    xcog = xcog_stacked(1);
    ycog = ycog_stacked(1);
    zcm_dist = zcm_dist_stacked(1);
    alpha = alpha_lone;
    beta = beta_lone;
    gamma = gamma_stacked(1);

    xcv = xcv_stacked(1);
    ycv = ycv_stacked(1);
    zcv = zcv_stacked(1);

    N_surf = N_surf_lone_cell{1};
    dA_surf = dA_surf_lone_cell{1};

    X_grid = X_grid_lone_cell{1};
    Y_grid = Y_grid_lone_cell{1};
    Z_grid = Z_grid_lone_cell{1};

    x_nodes = x_nodes_stacked(1,:);
    y_nodes = y_nodes_stacked(1,:);

    % rotate accordingly
    [x_proj_lone,y_proj_lone,H_lone,X_rot_lone,Y_rot_lone,Z_rot_lone,N_surf_rot_lone,dA_surf_rot_lone,X_top_lone,Y_top_lone,Z_top_lone,xcv_rot_lone,ycv_rot_lone,zcv_rot_lone] = RotateAndShift(X_grid,Y_grid,Z_grid,x_nodes,y_nodes,N_surf,dA_surf,alpha,beta,gamma,f,xcm,ycm,zcm_height_lone,to_plot_init,zcm_dist,xcog,ycog,thickness,is_flipped,xcv,ycv,zcv);

    X_rot_lone_cell{1} = X_rot_lone;
    Y_rot_lone_cell{1} = Y_rot_lone;
    Z_rot_lone_cell{1} = Z_rot_lone;
    X_top_lone_cell{1} = X_top_lone;
    Y_top_lone_cell{1} = Y_top_lone;
    Z_top_lone_cell{1} = Z_top_lone;
    N_surf_rot_lone_cell{1} = N_surf_rot_lone;
    dA_surf_rot_lone_cell{1} = dA_surf_rot_lone;

    % BEM solve for mu
    [mu_lone,rx_lone,ry_lone,nx_lone,ny_lone,w_lone,H_coll_lone,hs_lone] = mu_solve(x_proj_lone,y_proj_lone,N,H_lone,lc,gamma_E,1);

    % calc force and torque for the lone particle

    % calculate total force
    [Fz,Fx,Fy,Tx_buoyancy,Ty_buoyancy,Tz_buoyancy,...
    Fx_buoyancy,Fy_buoyancy,Fz_buoyancy,...
    Fx_meniscus, Fy_meniscus,Fz_meniscus] = CalcForce(X_rot_lone_cell{1},Y_rot_lone_cell{1},Z_rot_lone_cell{1},valid_mask_lone_cell{1},...
    N_surf_rot_lone_cell{1},dA_surf_rot_lone_cell{1},mu_lone,m,g,sigma,rho,w_lone,xcm_stacked(1),ycm_stacked(1),zcm_height_lone,...
    nx_lone,ny_lone,is_flipped,thickness,H_coll_lone,rx_lone,ry_lone,lc,hs_lone,0,0,0,X_top_lone_cell{1},...
    Y_top_lone_cell{1},Z_top_lone_cell{1},V_stacked(1),xcv_rot_lone,ycv_rot_lone);
    
    % calculate torques
    [Tx_men,Ty_men,Tz_men] = CalcTorque(mu_lone,w_lone,sigma,rx_lone,ry_lone,xcm_stacked(1),ycm_stacked(1),nx_lone,ny_lone,...
        H_coll_lone,zcm_height_lone,hs_lone,0);
    Tx = Tx_men + Tx_buoyancy;
    Ty = Ty_men + Ty_buoyancy;
    Tz = Tz_men + Tz_buoyancy;

    Fx_baseline = Fx;
    Fy_baseline = Fy;
    Fz_baseline = Fz;
    Tx_baseline = Tx;
    Ty_baseline = Ty;
    Tz_baseline = Tz;

    % compute Jacobian
    J_lone = Jacobian(1,x_nodes_stacked(1,:),y_nodes_stacked(1,:),xcm_stacked(1),ycm_stacked(1),...
        xcog_stacked(1,:),ycog_stacked(1,:),zcm_dist_stacked(1),gamma_stacked(1),m,zcm_height_lone,...
        N_surf_lone_cell,dA_surf_lone_cell,X_grid_lone_cell,Y_grid_lone_cell,Z_grid_lone_cell,...
        alpha_lone,beta_lone,x_proj_lone,y_proj_lone,H_lone,...
        X_rot_lone_cell,Y_rot_lone_cell,Z_rot_lone_cell,N_surf_rot_lone_cell,...
        dA_surf_rot_lone_cell,dz,dbeta,dalpha,f,to_plot_init,is_flipped,thickness,...
        N,lc,gamma_E,valid_mask_lone_cell,g,sigma,rho,0,0,0,0,xcv_stacked(1),ycv_stacked(1),zcv_stacked(1),...
        X_top_lone_cell,Y_top_lone_cell,Z_top_lone_cell,...
        xcv_rot_lone,ycv_rot_lone,zcv_rot_lone,V_stacked(1));


    iters = iters + 1;

    if max(abs(Fz_baseline)) <= tol && max(abs(Tx_baseline)) <= tol && max(abs(Ty_baseline)) <= tol
        break;
    end
    if iters >= max_iters
        break;
    end
    waitbar(iters/max_iters,hWaitbar,'Finding lone particle equilibrium...');

end
    

    
close(hWaitbar);
fprintf('\n INITIAL SINGLE-PARTICLE SOLVE FOR BASELINE\n');
fprintf(['\n *** Number of iterations: ',num2str(iters),' ***\n\n']);

fprintf('*** Baseline forces and torques for lone particle: ***\n\n');

fprintf(['Baseline force (z, vertical): ',num2str(Fz_baseline),'\n']);
fprintf(['Baseline torque (x, roll): ',num2str(Tx_baseline),'\n']);
fprintf(['Baseline torque (y, pitch): ',num2str(Ty_baseline),'\n\n']);
fprintf(['Baseline force (x, lateral): ',num2str(Fx_baseline),'\n']);
fprintf(['Baseline force (y, lateral): ',num2str(Fy_baseline),'\n']);
fprintf(['Baseline torque (z, yaw): ',num2str(Tz_baseline),'\n\n']);




if max(abs(mu_lone)) > 1
    disp('!!!!!!!!!!!!!!!!!!!!');
    fprintf('Small-slope approximation violated; results may be unphysical. Particle likely too dense, buoyant, inhomogeneous, or wavy\n');
    fprintf('Max(abs(mu)) = %.2f, should be < 1\n',max(abs(mu_lone)));
    disp('!!!!!!!!!!!!!!!!!!!!');
    disp('Press any key to continue anyway...');
    fprintf('\n');
    pause;
end


%% test BEM against manufactured solution

if check_manufactured

    % Choose a source point INSIDE the particle
    xs = xcm_stacked(1);
    ys = ycm_stacked(1);
    
    % Exact manufactured boundary condition at the BEM nodes
    
    dx = x_proj_lone' - xs;
    dy = y_proj_lone' - ys;
    r = sqrt(dx.^2 + dy.^2);
    
    H_test = besselk(0,r/lc);        

    [mu_BEM,rx_man,ry_man,nx_man,ny_man,~,~,~] = mu_solve(x_proj_lone,y_proj_lone,N,H_test,lc,gamma_E,1);

    mu_BEM = mu_BEM';
    
    % Exact mu at the collocation points
    
    dx = rx_man - xs;
    dy = ry_man - ys;
    r = sqrt(dx.^2 + dy.^2);
    
    mu_exact = -(1/lc) .* besselk(1,r/lc) .* ...
        (dx.*nx_man + dy.*ny_man) ./ r;
    
    % Error
    
    rel_error = norm(mu_BEM - mu_exact) / norm(mu_exact);
    
    fprintf('\n*** Check BEM solve against manufactured solution ***\n')
    fprintf('relative mu error = %.6e\n\n',rel_error);
end



%% now full assembly

% BEM solve for all particles
[mu_stacked,rx_stacked,ry_stacked,nx_stacked,ny_stacked,w_stacked,H_coll,hs_stacked] = mu_solve(x_proj_stacked,y_proj_stacked,N,H_stacked,lc,gamma_E,M);


Fz_stacked = zeros(M,1);
Fx_stacked = zeros(M,1);
Fy_stacked = zeros(M,1);
Tz_stacked = zeros(M,1);
Tx_stacked = zeros(M,1);
Ty_stacked = zeros(M,1);

fprintf('\n FINAL FULL ASSEMBLY SOLVE\n');

fprintf('*** Initial forces and torques: ***\n\n');

for p=1:M
    mu = mu_stacked((p-1)*N+1:p*N);
    m = m_stacked(p);
    X_rot = X_rot_stacked{p};
    Y_rot = Y_rot_stacked{p};
    Z_rot = Z_rot_stacked{p};
    X_top = X_top_stacked{p};
    Y_top = Y_top_stacked{p};
    Z_top = Z_top_stacked{p};
    xcv_rot = xcv_rot_stacked(p);
    ycv_rot = ycv_rot_stacked(p);
    zcv_rot = zcv_rot_stacked(p);
    N_surf_rot = N_surf_rot_stacked{p};
    dA_surf_rot = dA_surf_rot_stacked{p};
    w = w_stacked(p,:);
    xcm = xcm_stacked(p);
    ycm = ycm_stacked(p);
    zcm_height = zcm_height_stacked(p); % actual height of zcm
    r_x = rx_stacked(p,:);
    r_y = ry_stacked(p,:);
    nx = nx_stacked(p,:);
    ny = ny_stacked(p,:);
    H_coll_p = H_coll(p,:);
    x_proj = x_proj_stacked(p,:);
    y_proj = y_proj_stacked(p,:);
    hs = hs_stacked(p,:);
    gamma_rel = gamma_rel_stacked(p);
    V_particle = V_stacked(p);

    % calculate total force
    [Fz,Fx,Fy,Tx_buoyancy,Ty_buoyancy,Tz_buoyancy,...
    Fx_buoyancy,Fy_buoyancy,Fz_buoyancy,...
    Fx_meniscus, Fy_meniscus,Fz_meniscus] = CalcForce(X_rot,Y_rot,Z_rot,valid_mask_stacked{p},...
    N_surf_rot,dA_surf_rot,mu,m,g,sigma,rho,w,xcm,ycm,zcm_height,nx,ny,is_flipped,thickness,H_coll_p,r_x,r_y,lc,...
    hs,Fx_baseline,Fy_baseline,gamma_rel,X_top,Y_top,Z_top,V_particle,xcv_rot,ycv_rot);

    % calculate torques
    [Tx_men,Ty_men,Tz_men] = CalcTorque(mu,w,sigma,r_x,r_y,xcm,ycm,nx,ny,H_coll_p,zcm_height,hs,Tz_baseline);
    Tx = Tx_men + Tx_buoyancy;
    Ty = Ty_men + Ty_buoyancy;
    Tz = Tz_men + Tz_buoyancy;
    Fx_stacked(p) = Fx;
    Fy_stacked(p) = Fy;
    Fz_stacked(p) = Fz;
    Tx_stacked(p) = Tx;
    Ty_stacked(p) = Ty;
    Tz_stacked(p) = Tz;


    fprintf(['Particle ',num2str(p),':\n\n']);
    fprintf(['Initial force (z, vertical): ',num2str(Fz),'\n']);
    fprintf(['Initial torque (x, roll): ',num2str(Tx),'\n']);
    fprintf(['Initial torque (y, pitch): ',num2str(Ty),'\n\n']);
    fprintf(['Initial force (x, lateral): ',num2str(Fx),'\n']);
    fprintf(['Initial force (y, lateral): ',num2str(Fy),'\n']);
    fprintf(['Initial torque (z, yaw): ',num2str(Tz),'\n\n']);


end



%% see initial non-equilibrium meniscus solve

if view_init

    EvaluateDomain(1,min_x,min_y,max_x,max_y,grid_range,N_domain,x_proj_stacked,y_proj_stacked,...
        H_stacked,N,M,nx_stacked,ny_stacked,lc,w_stacked,rx_stacked,...
        ry_stacked,mu_stacked,X_rot_stacked,Y_rot_stacked,Z_rot_stacked,...
        xcm_stacked,ycm_stacked,zcm_height_stacked,limits,plotCM,...
        X_top_stacked, Y_top_stacked, Z_top_stacked, ...
        Fx_stacked, Fy_stacked, Fz_stacked, ...
        Tx_stacked, Ty_stacked, Tz_stacked,tol,L,H_coll,valid_mask_stacked,tol,tol,tol);

end



%% iteration to find equilibrium

Jacobians = cell(max_iters,1);

hWaitbar = waitbar(0,'Finding full assembly equilibrium...', 'Name', 'Status');


% compute Jacobian
J = Jacobian(M,x_nodes_stacked,y_nodes_stacked,xcm_stacked,ycm_stacked,...
    xcog_stacked,ycog_stacked,zcm_dist_stacked,gamma_stacked,m_stacked,zcm_height_stacked,...
    N_surf_stacked,dA_surf_stacked,X_grid_stacked,Y_grid_stacked,Z_grid_stacked,...
    alpha_stacked,beta_stacked,x_proj_stacked,y_proj_stacked,H_stacked,...
    X_rot_stacked,Y_rot_stacked,Z_rot_stacked,N_surf_rot_stacked,...
    dA_surf_rot_stacked,dz,dbeta,dalpha,f,to_plot_init,is_flipped,thickness,...
    N,lc,gamma_E,valid_mask_stacked,g,sigma,rho,Fx_baseline,Fy_baseline,Tz_baseline,gamma_rel_stacked,...
    xcv_stacked,ycv_stacked,zcv_stacked,X_top_stacked,Y_top_stacked,Z_top_stacked,...
    xcv_rot_stacked,ycv_rot_stacked,zcv_rot_stacked,V_stacked);


Jacobians{1} = J;


iters = 0;

while 1

    q = reshape([zcm_height_stacked, beta_stacked, alpha_stacked].', [], 1);
    R = reshape([Fz_stacked, Tx_stacked, Ty_stacked].', [], 1);

    dq = -J \ R;

    q = q + dq;

    zcm_height_stacked = q(1:3:end);
    beta_stacked = q(2:3:end);
    alpha_stacked = q(3:3:end);


    % shift and rotate particles accorindgly
    for p=1:M

        xcm = xcm_stacked(p);
        ycm = ycm_stacked(p);
        xcog = xcog_stacked(p);
        ycog = ycog_stacked(p);
        zcm_height = zcm_height_stacked(p);
        zcm_dist = zcm_dist_stacked(p);
        alpha = alpha_stacked(p);
        beta = beta_stacked(p);
        gamma = gamma_stacked(p);
        xcv = xcv_stacked(p);
        ycv = ycv_stacked(p);
        zcv = zcv_stacked(p);

        N_surf = N_surf_stacked{p};
        dA_surf = dA_surf_stacked{p};

        X_grid = X_grid_stacked{p};
        Y_grid = Y_grid_stacked{p};
        Z_grid = Z_grid_stacked{p};

        x_nodes = x_nodes_stacked(p,:);
        y_nodes = y_nodes_stacked(p,:);

        [x_proj,y_proj,H,X_rot,Y_rot,Z_rot,N_surf_rot,dA_surf_rot,X_top,Y_top,Z_top,xcv_rot,ycv_rot,zcv_rot] = RotateAndShift(X_grid,Y_grid,Z_grid,x_nodes,y_nodes,N_surf,dA_surf,alpha,beta,gamma,f,xcm,ycm,zcm_height,to_plot_init,zcm_dist,xcog,ycog,thickness,is_flipped,xcv,ycv,zcv);

        x_proj_stacked(p,:) = x_proj;
        y_proj_stacked(p,:) = y_proj;
        H_stacked(:,p) = H;
        X_rot_stacked{p} = X_rot;
        Y_rot_stacked{p} = Y_rot;
        Z_rot_stacked{p} = Z_rot;
        X_top_stacked{p} = X_top;
        Y_top_stacked{p} = Y_top;
        Z_top_stacked{p} = Z_top;
        xcv_rot_stacked(p) = xcv_rot;
        ycv_rot_stacked(p) = ycv_rot;
        zcv_rot_stacked(p) = zcv_rot;
        N_surf_rot_stacked{p} = N_surf_rot;
        dA_surf_rot_stacked{p} = dA_surf_rot;

    end

    % BEM solve for mu
    [mu_stacked,rx_stacked,ry_stacked,nx_stacked,ny_stacked,w_stacked,H_coll,hs_stacked] = mu_solve(x_proj_stacked,y_proj_stacked,N,H_stacked,lc,gamma_E,M);


    % calc force and torque for each particle
    for p=1:M
        mu = mu_stacked((p-1)*N+1:p*N);
        m = m_stacked(p);
        X_rot = X_rot_stacked{p};
        Y_rot = Y_rot_stacked{p};
        Z_rot = Z_rot_stacked{p};
        X_top = X_top_stacked{p};
        Y_top = Y_top_stacked{p};
        Z_top = Z_top_stacked{p};
        N_surf_rot = N_surf_rot_stacked{p};
        dA_surf_rot = dA_surf_rot_stacked{p};
        w = w_stacked(p,:);
        xcm = xcm_stacked(p);
        ycm = ycm_stacked(p);
        zcm_height = zcm_height_stacked(p);
        r_x = rx_stacked(p,:);
        r_y = ry_stacked(p,:);
        nx = nx_stacked(p,:);
        ny = ny_stacked(p,:);
        H = H_stacked(:,p);
        H_coll_p = H_coll(p,:);
        x_proj = x_proj_stacked(p,:);
        y_proj = y_proj_stacked(p,:);
        hs = hs_stacked(p,:);
        gamma_rel = gamma_rel_stacked(p);
        V_particle = V_stacked(p);
        xcv_rot = xcv_rot_stacked(p);
        ycv_rot = ycv_rot_stacked(p);


        % calculate total force
        [Fz,Fx,Fy,Tx_buoyancy,Ty_buoyancy,Tz_buoyancy,...
                Fx_buoyancy,Fy_buoyancy,Fz_buoyancy,...
                Fx_meniscus, Fy_meniscus,Fz_meniscus] = CalcForce(X_rot,Y_rot,Z_rot,valid_mask_stacked{p},N_surf_rot,...
                dA_surf_rot,mu,m,g,sigma,rho,w,xcm,ycm,zcm_height,nx,ny,is_flipped,thickness,H_coll_p,...
                r_x,r_y,lc,hs,Fx_baseline,Fy_baseline,gamma_rel,X_top,Y_top,Z_top,V_particle,xcv_rot,ycv_rot);


        % calculate torques
        [Tx_men,Ty_men,Tz_men] = CalcTorque(mu,w,sigma,r_x,r_y,xcm,ycm,nx,ny,H_coll_p,zcm_height,hs,Tz_baseline);
        Tx = Tx_men + Tx_buoyancy;
        Ty = Ty_men + Ty_buoyancy;
        Tz = Tz_men + Tz_buoyancy;
        Fx_stacked(p) = Fx;
        Fy_stacked(p) = Fy;
        Fz_stacked(p) = Fz;
        Tx_stacked(p) = Tx;
        Ty_stacked(p) = Ty;
        Tz_stacked(p) = Tz;

        X_top_stacked{p} = X_top;
        Y_top_stacked{p} = Y_top;
        Z_top_stacked{p} = Z_top;



    end


    % compute Jacobian
    J = Jacobian(M,x_nodes_stacked,y_nodes_stacked,xcm_stacked,ycm_stacked,...
    xcog_stacked,ycog_stacked,zcm_dist_stacked,gamma_stacked,m_stacked,zcm_height_stacked,...
    N_surf_stacked,dA_surf_stacked,X_grid_stacked,Y_grid_stacked,Z_grid_stacked,...
    alpha_stacked,beta_stacked,x_proj_stacked,y_proj_stacked,H_stacked,...
    X_rot_stacked,Y_rot_stacked,Z_rot_stacked,N_surf_rot_stacked,...
    dA_surf_rot_stacked,dz,dbeta,dalpha,f,to_plot_init,is_flipped,thickness,...
    N,lc,gamma_E,valid_mask_stacked,g,sigma,rho,Fx_baseline,Fy_baseline,Tz_baseline,gamma_rel_stacked,...
    xcv_stacked,ycv_stacked,zcv_stacked,X_top_stacked,Y_top_stacked,Z_top_stacked,...
    xcv_rot_stacked,ycv_rot_stacked,zcv_rot_stacked,V_stacked);

    iters = iters + 1;

    Jacobians{iters+1} = J;

    if max(abs(Fz_stacked)) <= tol && max(abs(Tx_stacked)) <= tol && max(abs(Ty_stacked)) <= tol
        break;
    end
    if iters >= max_iters
        break;
    end
    waitbar(iters/max_iters,hWaitbar,'Finding full assembly equilibrium...');

end
    

    
close(hWaitbar);
fprintf(['\n *** Number of iterations: ',num2str(iters),' ***\n\n']);

fprintf('*** Final forces and torques: ***\n\n');

for p=1:M
    Fx = Fx_stacked(p);
    Fy = Fy_stacked(p);
    Fz = Fz_stacked(p);
    Tx = Tx_stacked(p);
    Ty = Ty_stacked(p);
    Tz = Tz_stacked(p);
    fprintf(['Particle ',num2str(p),':\n\n']);
    fprintf(['Final force (z, vertical): ',num2str(Fz),'\n']);
    fprintf(['Final torque (x, roll): ',num2str(Tx),'\n']);
    fprintf(['Final torque (y, pitch): ',num2str(Ty),'\n\n']);
    fprintf(['Final force (x, lateral): ',num2str(Fx),'\n']);
    fprintf(['Final force (y, lateral): ',num2str(Fy),'\n']);
    fprintf(['Final torque (z, yaw): ',num2str(Tz),'\n\n']);
end

fprintf('*** Check for conservation of angular momentum (zero net torque) and zero net force: ***\n\n');
[residual_torque,residual_force_x,residual_force_y] = CheckGlobalTorque(xcm_stacked, ycm_stacked, Fx_stacked, Fy_stacked, Tz_stacked);

%% look at final equilibrium solution

EvaluateDomain(0,min_x,min_y,max_x,max_y,grid_range,N_domain,x_proj_stacked,y_proj_stacked,...
        H_stacked,N,M,nx_stacked,ny_stacked,lc,w_stacked,rx_stacked,...
        ry_stacked,mu_stacked,X_rot_stacked,Y_rot_stacked,Z_rot_stacked,...
        xcm_stacked,ycm_stacked,zcm_height_stacked,limits,plotCM,...
        X_top_stacked, Y_top_stacked, Z_top_stacked, ...
        Fx_stacked, Fy_stacked, Fz_stacked, ...
        Tx_stacked, Ty_stacked, Tz_stacked,tol,L,H_coll,valid_mask_stacked,residual_torque,residual_force_x,residual_force_y);


%% plot Jacobians

Jacobians = Jacobians(~cellfun('isempty', Jacobians));
Jacobians_scaled = Jacobians;


F0 = m*g;
Sq = diag(repmat([L, 1, 1], 1, M));
Sr = diag(repmat([F0, F0*L, F0*L], 1, M));

for i=1:length(Jacobians)
    Jacobians_scaled{i} = Sr \ Jacobians{i} * Sq;
end


PlotJacobian(M,Jacobians_scaled{1},'dimensionless initial ');

PlotJacobian(M,Jacobians_scaled{end},'dimensionless final ');


%% show Jacobian evolution over time

diagonal_metric = zeros(iters+1,1);
symmetry_metric = zeros(iters+1,1);

for i = 0:iters

    J = Jacobians_scaled{i+1};

    J_offdiag = J - diag(diag(J));

    % relative size of off-diagonal entries
    diagonal_metric(i+1) = norm(J_offdiag,'fro') / norm(J,'fro');

    % relative size of antisymmetric part
    symmetry_metric(i+1) = norm(J - J.','fro') / norm(J,'fro');

end

figure;

subplot(2,1,1);
plot(0:iters,diagonal_metric,'o-','LineWidth',1.5);
xlabel('Iteration','Interpreter','latex');
ylabel('$\|J-\rm{diag}(J)\|_F/\|J\|_F$', ...
    'Interpreter','latex');
grid on;
title('how diagonal?');

subplot(2,1,2);
plot(0:iters,symmetry_metric,'o-','LineWidth',1.5);
xlabel('Iteration','Interpreter','latex');
ylabel('$\|J-J^T\|_F/\|J\|_F$', ...
    'Interpreter','latex');
grid on;
title('how symmetric?');


%% look at eigenvalues & eigenvectors

[V,D] = eig(Jacobians_scaled{end});


figure;
imagesc(V);
axis image;

% Colormap
colorbar;
colormap(redwhiteblue(-1,1));

% Tick locations: one for each variable within each particle
ticks = 1:3*M;

labels_x = cell(3*M,1);
labels_y = cell(3*M,1);
for k = 1:3*M
    labels_x{k} = num2str(round(D(k,k),2));
end

for k=1:M
    labels_y{3*k-2} = sprintf('$\\delta z^{%d}$',k);
    labels_y{3*k-1} = sprintf('$\\delta\\beta^{%d}$',k);
    labels_y{3*k}   = sprintf('$\\delta\\alpha^{%d}$',k);
end

xticks(ticks);
yticks(ticks);
xticklabels(labels_x);
yticklabels(labels_y);

set(gca, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',12);

xlabel('$\lambda$','Interpreter','latex','FontSize',16);
ylabel('$\mathbf{V}$','Interpreter','latex','FontSize',16);

title('eigenvectors, $\mathbf{V}$', ...
    'Interpreter','latex','FontSize',16);

hold on;

for k = 3:3:3*M-1
    yline(k + 0.5,'k-','LineWidth',1);
end

for k=M:M:3*M-1
    xline(k + 0.5,'k-','LineWidth',1);
end

clim([-1 1]);

hold off;