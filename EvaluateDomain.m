function EvaluateDomain(is_init,min_x,min_y,max_x,max_y,grid_range,N_domain,x_proj_stacked,...
    y_proj_stacked,H_stacked,N,M,nx_stacked,ny_stacked,lc,w_stacked,...
    rx_stacked,ry_stacked,mu_stacked,X_rot_stacked,Y_rot_stacked,...
    Z_rot_stacked,xcm_stacked,ycm_stacked,zcm_height_stacked,limits,plotCM,...
    X_top_stacked, Y_top_stacked, Z_top_stacked, ...
    Fx_stacked, Fy_stacked, Fz_stacked,...
    Tx_stacked, Ty_stacked, Tz_stacked,tol,L,H_coll,valid_mask_stacked,residual_torque,residual_force_x,residual_force_y)

center = [(max_x+min_x)/2, (max_y+min_y)/2];

[X, Y] = meshgrid(linspace(center(1)-grid_range/2, center(1)+grid_range/2, N_domain),linspace(center(2)-grid_range/2, center(2)+grid_range/2, N_domain));
Z = NaN(size(X)); 

grid_spacing = X(1,2) - X(1,1);

if ~is_init
    hWaitbar = waitbar(0,'Evaluating final fluid domain...', 'Name', 'Status');
else
    hWaitbar = waitbar(0,'Evaluating initial fluid domain...', 'Name', 'Status');
end

for q = 1:numel(X)
    next_q = false;
    h_q = 0;
    % skip points inside all bodies
    for k=1:M
        x_proj = x_proj_stacked(k,:);
        y_proj = y_proj_stacked(k,:);
        r_x = rx_stacked(k,:);
        r_y = ry_stacked(k,:);
        nx = nx_stacked(k,:);
        ny = ny_stacked(k,:);
        H = H_stacked(:,k);
        H_coll_p = H_coll(k,:);
        w = w_stacked(k,:);
        mu = mu_stacked((k-1)*N+1:k*N);

        x_c = mean(x_proj);
        y_c = mean(y_proj);
        
        dx = x_proj - x_c;
        dy = y_proj - y_c;
        
        r = sqrt(dx.^2 + dy.^2);
        
        x_proj_expanded = x_c + (r + grid_spacing).*dx./r;
        y_proj_expanded = y_c + (r + grid_spacing).*dy./r;

        % if it's really close to or inside of body k
        if inpolygon(X(q), Y(q), x_proj_expanded, y_proj_expanded)
            if inpolygon(X(q), Y(q), x_proj, y_proj)
                next_q = true;
                break;
            else
                rq = [X(q); Y(q)];
                r = [r_x; r_y];
                dist = vecnorm(r-rq,2,1);
                [~,i_near] = min(dist);
                %Z(q) = H(i_near);
                Z(q) = H_coll_p(i_near);
                next_q = true;
                break;
            end
        end
        if next_q
            continue;
        end
        
        % otherwise, full evaluation
        for i = 1:N
            dist_x = r_x(i) - X(q);
            dist_y = r_y(i) - Y(q);
            r = sqrt(dist_x^2 + dist_y^2);
            
            ndotr = (dist_x * nx(i) + dist_y * ny(i)) / r;
            
            % q is the evaluation point, i is the integration boundary point
            B_iq = (1/(2*pi*lc)) * besselk(1, r/lc) * ndotr;
            G_iq = -(1/(2*pi)) * besselk(0, r/lc);
            
            % h_k = sum(H_j * B_jk * w_j - mu_j * G_jk * w_j)
            h_q = h_q + (B_iq * H_coll_p(i) - G_iq * mu(i)) * w(i);
        end
    end
    if ~next_q
        Z(q) = h_q;
    end
    if ~is_init
        waitbar(q/numel(X),hWaitbar,'Evaluating final fluid domain...');
    else
        waitbar(q/numel(X),hWaitbar,'Evaluating initial fluid domain...');
    end
end

close(hWaitbar);

%% Visualization
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
colormap(jet); 
hold on;


for k=1:M
    surf(X_rot_stacked{k}, Y_rot_stacked{k}, Z_rot_stacked{k}, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    hold on
    if plotCM
        scatter3(xcm_stacked(k),ycm_stacked(k),zcm_height_stacked(k),200,'filled','red');
    end
    scatter3(x_proj_stacked(k,:),y_proj_stacked(k,:),H_stacked(:,k),'MarkerEdgeColor','black');
    text(xcm_stacked(k), ycm_stacked(k), zcm_height_stacked(k),num2str(k)); 
end



% Formatting
if ~is_init
    title('Final meniscus solution', 'FontSize', 14);
else
    title('Initial meniscus solution', 'FontSize', 14);
end
xlabel('x (cm)', 'FontSize', 14); ylabel('y (cm)', 'FontSize', 14); zlabel('h (cm)', 'FontSize', 14);
axis equal; 
grid on;
view(3)
box on

clim(limits);


colorbar


hLegendF = gobjects(0);
hLegendT = gobjects(0);



for p=1:M
    Xc = X_rot_stacked{p};
    Yc = Y_rot_stacked{p};
    Zc = Z_rot_stacked{p};

    
    Xt = X_top_stacked{p};
    Yt = Y_top_stacked{p};
    Zt = Z_top_stacked{p};
    
    valid = ~isnan(Zc) & ~isnan(Zt);
    
    % Find the boundary of the valid region
    B = bwboundaries(valid, 'noholes');
    
    for ib = 1:numel(B)
        b = B{ib};
    
        ii = b(:,1);
        jj = b(:,2);
    
        % Bottom boundary
        Xb = Xc(sub2ind(size(Xc),ii,jj));
        Yb = Yc(sub2ind(size(Yc),ii,jj));
        Zb = Zc(sub2ind(size(Zc),ii,jj));
    
        % Top boundary
        Xt_b = Xt(sub2ind(size(Xt),ii,jj));
        Yt_b = Yt(sub2ind(size(Yt),ii,jj));
        Zt_b = Zt(sub2ind(size(Zt),ii,jj));
    
        % Connect the two corresponding boundaries
        surf([Xb, Xt_b], ...
             [Yb, Yt_b], ...
             [Zb, Zt_b], ...
             'FaceColor', 'k', ...
             'EdgeColor', 'none', ...
             'FaceAlpha', 0.4);
    end
    % Top boundary only
    surf(Xt, Yt, Zt, ...
    'FaceColor', 'k', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.4);

    xcm = xcm_stacked(p);
    ycm = ycm_stacked(p);
    Fx = Fx_stacked(p);
    Fy = Fy_stacked(p);
    Fz = Fz_stacked(p);
    Tx = Tx_stacked(p);
    Ty = Ty_stacked(p);
    Tz = Tz_stacked(p);
    zcm_height = zcm_height_stacked(p);

    % only plot force/torque vectors if they're decently bigger than the residual system force/torque

    Fmag = sqrt(Fx^2 + Fy^2 + Fz^2);

    if floor(log10(abs(Fmag))) > max(floor(log10(abs(tol))),max(floor(log10(abs(residual_force_x))),floor(log10(abs(residual_force_y)))))
        % Force vector
        hF = quiver3(xcm, ycm, zcm_height, L*Fx/Fmag, L*Fy/Fmag, L*Fz/Fmag, ...
            0, 'LineWidth', 4, 'MaxHeadSize', 1,'DisplayName',['|F|=',num2str(Fmag),' dynes']);
        hLegendF(end+1) = plot3(NaN, NaN, NaN, ...
            'Color', hF.Color, 'LineWidth', 3, ...
            'DisplayName', hF.DisplayName);
    end


    Tmag = sqrt(Tx^2 + Ty^2 + Tz^2);

    if floor(log10(abs(Tmag))) > max(floor(log10(abs(tol))),floor(log10(abs(residual_torque))))
        % Torque vector
        hT = quiver3(xcm, ycm, zcm_height, L*Tx/Tmag, L*Ty/Tmag, L*Tz/Tmag, ...
            0, 'LineWidth', 4, 'MaxHeadSize', 1,'DisplayName',['|T|=',num2str(Tmag),' dynes x cm']);
        hLegendT(end+1) = plot3(NaN, NaN, NaN, ...
            'Color', hT.Color, 'LineWidth', 3, ...
            'DisplayName', hT.DisplayName);
    end
    
    axis equal
     



end


if ~isempty(hLegendF) || ~isempty(hLegendT)
    legend([hLegendF, hLegendT], ...
        'Location', 'best', ...
        'FontSize', 14);
else
    title('Final meniscus solution. Negligible force & torque');
end



end

