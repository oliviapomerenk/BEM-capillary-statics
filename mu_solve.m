function [mu_stacked,rx_stacked,ry_stacked,nx_stacked,ny_stacked,w_stacked,H_coll,hs_stacked] = mu_solve(x_nodes_stacked,y_nodes_stacked,N,H_stacked,lc,gamma_E,M)

% build weights and normals
w_stacked = zeros(M,N);
nx_stacked = zeros(M,N);
ny_stacked = zeros(M,N);
rx_stacked = zeros(M,N);
ry_stacked = zeros(M,N);



for k=1:M

    x_nodes = x_nodes_stacked(k,:);
    y_nodes = y_nodes_stacked(k,:); 

    for i = 1:N
        x1 = x_nodes(i); y1 = y_nodes(i);
        % periodicity
        if i == N
            x2 = x_nodes(1); y2 = y_nodes(1);
        else
            x2 = x_nodes(i+1); y2 = y_nodes(i+1);
        end

        % Collocation points r_i (midpoint of each segment)
        rx_stacked(k,i) = (x1 + x2) / 2;
        ry_stacked(k,i) = (y1 + y2) / 2;

        dx = x2 - x1; 
        dy = y2 - y1;
        w_stacked(k,i) = sqrt(dx^2 + dy^2); % weight is length of segment

        % Normal vector pointing INWARD to the object (OUTWARD from fluid)
        nx_stacked(k,i) = -dy / w_stacked(k,i);
        ny_stacked(k,i) = dx / w_stacked(k,i);
    end

end


H_coll = zeros(M,N);

hs_stacked = zeros(M,N); % Initialize the tangential slope array

for k = 1:M

    H_nodes = H_stacked(:,k);

    for i = 1:N

        if i == N
            % Periodic last segment
            H_coll(k,i) = 0.5 * (H_nodes(i) + H_nodes(1));
            dh = H_nodes(1) - H_nodes(i);
        else
            H_coll(k,i) = 0.5 * (H_nodes(i) + H_nodes(i+1));
            % dh/ds for standard segments
            dh = H_nodes(i+1) - H_nodes(i);
        end

        % The tangential slope is exactly dh divided by the segment arc length
        hs_stacked(k,i) = dh / w_stacked(k,i);

    end
end


A4D = zeros(M,M,N,N);
B4D = zeros(M,M,N,N);

for p=1:M
    B_p = zeros(N,N);
    G_p = zeros(N,N);
    for i = 1:N
        for j = 1:N
            if i == j
                % diagonal
                B_p(j,i) = 0;
                %G_p(j,i) = (1/(2*pi)) * w_stacked(p,i) * (log(w_stacked(p,i)/(4*lc)) + gamma_E - 1);
                G_p(j,i) = (1/(2*pi)) * (log(w_stacked(p,i)/(4*lc)) + gamma_E - 1);
            else
                % integration
                dist_x = rx_stacked(p,j) - rx_stacked(p,i);
                dist_y = ry_stacked(p,j) - ry_stacked(p,i);
                r = sqrt(dist_x^2 + dist_y^2);

                % (r_j - r_i) / |r_j - r_i| dot n(r_j)
                ndotr = (dist_x * nx_stacked(p,j) + dist_y * ny_stacked(p,j)) / r;

                B_p(j,i) = (1/(2*pi*lc)) * besselk(1, r/lc) * ndotr;
                G_p(j,i) = -(1/(2*pi)) * besselk(0, r/lc);
            end
            A4D(p,p,i,j) = G_p(j,i) * w_stacked(p,j);
        end
    end
    B4D(p,p,:,:) = B_p;

    for k=1:M
        if k ~= p
            B_kp = zeros(N,N);
            G_kp = zeros(N,N);
            for i=1:N
                for j=1:N
                    % integration
                    dist_x = rx_stacked(k,j) - rx_stacked(p,i);
                    dist_y = ry_stacked(k,j) - ry_stacked(p,i);
                    r = sqrt(dist_x^2 + dist_y^2);

                    % (r_j - r_i) / |r_j - r_i| dot n(r_j)
                    ndotr = (dist_x * nx_stacked(k,j) + dist_y * ny_stacked(k,j)) / r;

                    B_kp(j,i) = (1/(2*pi*lc)) * besselk(1, r/lc) * ndotr;
                    G_kp(j,i) = -(1/(2*pi)) * besselk(0, r/lc);
                    % Construct the A_pk block matrix
                    A4D(p,k,i,j) = G_kp(j,i) * w_stacked(k,j);
                end
            end
            B4D(k,p,:,:) = B_kp;
        end
    end
end



% build A
A = zeros(M*N, M*N);
for p = 1:M
    rows = (p-1)*N + (1:N);
    for k = 1:M
        cols = (k-1)*N + (1:N);
        A(rows,cols) = squeeze(A4D(p,k,:,:));
    end
end

% build b

b = zeros(M*N,1);

for p = 1:M
    rows = (p-1)*N + (1:N);

    bp = -0.5 * H_coll(p,:).';

    for k = 1:M
        Bkp = squeeze(B4D(k,p,:,:));  % B_ji^(k,p)

        Hkwk = H_coll(k,:).' .* w_stacked(k,:).';

        bp = bp + Bkp.' * Hkwk;
    end

    b(rows) = bp;
end


mu_stacked = A \ b;

end