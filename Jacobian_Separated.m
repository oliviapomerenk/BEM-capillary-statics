function [J, J_cap, J_hydro] = Jacobian_Separated(M,x_nodes_stacked,y_nodes_stacked,xcm_stacked,ycm_stacked,...
    xcog_stacked,ycog_stacked,zcm_dist_stacked,gamma_stacked,m_stacked,zcm_height_stacked,...
    N_surf_stacked,dA_surf_stacked,X_grid_stacked,Y_grid_stacked,Z_grid_stacked,...
    alpha_stacked,beta_stacked,x_proj_stacked,y_proj_stacked,H_stacked,...
    X_rot_stacked,Y_rot_stacked,Z_rot_stacked,N_surf_rot_stacked,...
    dA_surf_rot_stacked,dz,dbeta,dalpha,f,to_plot_init,is_flipped,thickness,...
    N,lc,gamma_E,valid_mask_stacked,g,sigma,rho,Fx_baseline,Fy_baseline,Tz_baseline,...
    gamma_rel_stacked,xcv_stacked,ycv_stacked,zcv_stacked,X_top_stacked,Y_top_stacked,Z_top_stacked,...
    xcv_rot_stacked,ycv_rot_stacked,zcv_rot_stacked,V_stacked)



J_cap_cell = cell(M,M);
J_cap_plus = zeros(3,3,M,M); % plus block in Jacobian matrix
J_cap_minus = zeros(3,3,M,M); % minus block in Jacobian matrix

J_hydro_cell = cell(M,M);
J_hydro_plus = zeros(3,3,M,M); % plus block in Jacobian matrix
J_hydro_minus = zeros(3,3,M,M); % minus block in Jacobian matrix


for k=1:M

    % for each perturbation, rotate & shift particle k accordingly -> full BEM solve -> calc force on all particles

    for i=1:6

        % reset quantities
        zcm_height = zcm_height_stacked(k);
        alpha = alpha_stacked(k);
        beta = beta_stacked(k);

        % everything for particle k
        x_nodes = x_nodes_stacked(k,:);
        y_nodes = y_nodes_stacked(k,:);
        xcm = xcm_stacked(k);
        ycm = ycm_stacked(k);
        xcog = xcog_stacked(k);
        ycog = ycog_stacked(k);
        zcm_dist = zcm_dist_stacked(k);
        gamma = gamma_stacked(k);
        m = m_stacked(k);
        gamma_rel = gamma_rel_stacked(k);
        xcv = xcv_stacked(k);
        ycv = ycv_stacked(k);
        zcv = zcv_stacked(k);
    
        % normals and area elements for particle k
        N_surf = N_surf_stacked{k};
        dA_surf = dA_surf_stacked{k};
    
        % baseline grid for particle k
        X_grid = X_grid_stacked{k};
        Y_grid = Y_grid_stacked{k};
        Z_grid = Z_grid_stacked{k};
        
        % make copies so we can perturb them without overwriting anything
        x_proj_pert = x_proj_stacked;
        y_proj_pert = y_proj_stacked;
        H_pert = H_stacked;
        X_rot_pert = X_rot_stacked;
        Y_rot_pert = Y_rot_stacked;
        Z_rot_pert = Z_rot_stacked;
        X_top_pert = X_top_stacked;
        Y_top_pert = Y_top_stacked;
        Z_top_pert = Z_top_stacked;
        N_surf_rot_pert = N_surf_rot_stacked;
        dA_surf_rot_pert = dA_surf_rot_stacked;
        zcm_height_pert = zcm_height_stacked;
        xcv_rot_pert = xcv_rot_stacked;
        ycv_rot_pert = ycv_rot_stacked;

        % perturb
        if i==1
            zcm_height = zcm_height + dz;
            zcm_height_pert(k) = zcm_height;
        elseif i==2
            zcm_height = zcm_height - dz;
            zcm_height_pert(k) = zcm_height;
        elseif i==3
            beta = beta + dbeta;
        elseif i==4
            beta = beta - dbeta;
        elseif i==5
            alpha = alpha + dalpha;
        else
            alpha = alpha - dalpha;
        end

        % rotate accordingly
        [x_proj,y_proj,H,X_rot,Y_rot,Z_rot,N_surf_rot,dA_surf_rot,X_top,Y_top,Z_top,xcv_rot,ycv_rot,zcv_rot] = RotateAndShift(X_grid,Y_grid,Z_grid,x_nodes,y_nodes,N_surf,dA_surf,alpha,beta,gamma,f,xcm,ycm,zcm_height,to_plot_init,zcm_dist,xcog,ycog,thickness,is_flipped,xcv,ycv,zcv);
        
        x_proj_pert(k,:) = x_proj;
        y_proj_pert(k,:) = y_proj;
        H_pert(:,k) = H;
        X_rot_pert{k} = X_rot;
        Y_rot_pert{k} = Y_rot;
        Z_rot_pert{k} = Z_rot;
        X_top_pert{k} = X_top;
        Y_top_pert{k} = Y_top;
        Z_top_pert{k} = Z_top;
        N_surf_rot_pert{k} = N_surf_rot;
        dA_surf_rot_pert{k} = dA_surf_rot;
        xcv_rot_pert(k) = xcv_rot;
        ycv_rot_pert(k) = ycv_rot;

        % BEM solve for perturbed system
        [mu_stacked,rx_stacked,ry_stacked,nx_stacked,ny_stacked,w_stacked,H_coll,hs_stacked] = mu_solve(x_proj_pert,y_proj_pert,N,H_pert,lc,gamma_E,M);
    

        % calc force on all particles based on this particular perturbation
        for p=1:M

            mu = mu_stacked((p-1)*N+1:p*N);
            w = w_stacked(p,:);
            r_x = rx_stacked(p,:);
            r_y = ry_stacked(p,:);
            nx = nx_stacked(p,:);
            ny = ny_stacked(p,:);
            H_coll_p = H_coll(p,:);
            X_rot = X_rot_pert{p};
            Y_rot = Y_rot_pert{p};
            Z_rot = Z_rot_pert{p};
            X_top = X_top_pert{p};
            Y_top = Y_top_pert{p};
            Z_top = Z_top_pert{p};
            N_surf_rot = N_surf_rot_pert{p};
            dA_surf_rot = dA_surf_rot_pert{p};
            xcm = xcm_stacked(p);
            ycm = ycm_stacked(p);
            zcm_height = zcm_height_pert(p);
            hs = hs_stacked(p,:);
            xcv_rot = xcv_rot_pert(p);
            ycv_rot = ycv_rot_pert(p);
            V_particle = V_stacked(p);

    
            % calculate total force
            [Fz,Fx,Fy,Tx_buoyancy,Ty_buoyancy,Tz_buoyancy,...
            Fx_buoyancy,Fy_buoyancy,Fz_buoyancy,...
            Fx_meniscus, Fy_meniscus,Fz_meniscus] = CalcForce(X_rot,Y_rot,Z_rot,valid_mask_stacked{p},N_surf_rot,dA_surf_rot,mu,m,g,sigma,...
            rho,w,xcm,ycm,zcm_height,nx,ny,is_flipped,thickness,H_coll_p,r_x,r_y,lc,hs,Fx_baseline,Fy_baseline,gamma_rel,...
            X_top,Y_top,Z_top,V_particle,xcv_rot,ycv_rot);
            
            % calculate torques
            [Tx_men,Ty_men,Tz_men] = CalcTorque(mu,w,sigma,r_x,r_y,xcm,ycm,nx,ny,H_coll_p,zcm_height,hs,Tz_baseline);
            Tx = Tx_men + Tx_buoyancy;
            Ty = Ty_men + Ty_buoyancy;
            Tz = Tz_men + Tz_buoyancy;

            if i==1 
                J_cap_plus(:,1,p,k) = [Fz_meniscus; Tx_men; Ty_men];
                J_hydro_plus(:,1,p,k) = [Fz_buoyancy; Tx_buoyancy; Ty_buoyancy];
            elseif i==2
                J_cap_minus(:,1,p,k) = [Fz_meniscus; Tx_men; Ty_men];
                J_hydro_minus(:,1,p,k) = [Fz_buoyancy; Tx_buoyancy; Ty_buoyancy];
            elseif i==3
                J_cap_plus(:,2,p,k) = [Fz_meniscus; Tx_men; Ty_men];
                J_hydro_plus(:,2,p,k) = [Fz_buoyancy; Tx_buoyancy; Ty_buoyancy];
            elseif i==4
                J_cap_minus(:,2,p,k) = [Fz_meniscus; Tx_men; Ty_men];
                J_hydro_minus(:,2,p,k) = [Fz_buoyancy; Tx_buoyancy; Ty_buoyancy];
            elseif i==5
                J_cap_plus(:,3,p,k) = [Fz_meniscus; Tx_men; Ty_men];
                J_hydro_plus(:,3,p,k) = [Fz_buoyancy; Tx_buoyancy; Ty_buoyancy];
            else
                J_cap_minus(:,3,p,k) = [Fz_meniscus; Tx_men; Ty_men];
                J_hydro_minus(:,3,p,k) = [Fz_buoyancy; Tx_buoyancy; Ty_buoyancy];
            end
        end
    end

end

for k=1:M
    for p=1:M
        Jpk_cap = zeros(3,3);
        Jpk_hydro = zeros(3,3);

        Jpk_cap_plus = J_cap_plus(:,:,p,k);
        Jpk_cap_minus = J_cap_minus(:,:,p,k);
        Jpk_hydro_plus = J_hydro_plus(:,:,p,k);
        Jpk_hydro_minus = J_hydro_minus(:,:,p,k);

        Jpk_cap(:,1) = (Jpk_cap_plus(:,1) - Jpk_cap_minus(:,1)) / (2*dz);
        Jpk_cap(:,2) = (Jpk_cap_plus(:,2) - Jpk_cap_minus(:,2)) / (2*dbeta);
        Jpk_cap(:,3) = (Jpk_cap_plus(:,3) - Jpk_cap_minus(:,3)) / (2*dalpha);

        Jpk_hydro(:,1) = (Jpk_hydro_plus(:,1) - Jpk_hydro_minus(:,1)) / (2*dz);
        Jpk_hydro(:,2) = (Jpk_hydro_plus(:,2) - Jpk_hydro_minus(:,2)) / (2*dbeta);
        Jpk_hydro(:,3) = (Jpk_hydro_plus(:,3) - Jpk_hydro_minus(:,3)) / (2*dalpha);

        J_cap_cell{p,k} = Jpk_cap;
        J_hydro_cell{p,k} = Jpk_hydro;
    end
end

J_cap = cell2mat(J_cap_cell); % full 3M x 3M Jacobian
J_hydro = cell2mat(J_hydro_cell); % full 3M x 3M Jacobian

J = J_cap + J_hydro;


end