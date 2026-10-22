function PlotJacobian(M,J_initial,J_cap_initial,J_hydro_initial,...
                         J_final,J_cap_final,J_hydro_final,title_string)

% used ChatGPT for colorbar stuff

figure;

% Matrices and titles
J_all = {J_initial, J_cap_initial, J_hydro_initial; ...
         J_final,   J_cap_final,   J_hydro_final};

J_titles = { ...
    '$\mathbf{\tilde{J}}=\mathbf{\tilde{J}}_{\rm cap}+\mathbf{\tilde{J}}_{\rm hydro}$', ...
    '$\mathbf{\tilde{J}}_{\rm cap}$', ...
    '$\mathbf{\tilde{J}}_{\rm hydro}$'};

% Symmetric color limits, always centered at zero
clim_max = max(abs(J_initial(:)));
color_limits = [-clim_max, clim_max];

% Colormap: dark red -> white -> light blue -> light purple
n = 256;
nneg = n/2;
npos = n - nneg;

cmap_hot = hot(256);
cmap_neg = interp1(linspace(0,1,256), cmap_hot, ...
                   linspace(0.3,1,nneg));

white  = [1.00 1.00 1.00];
blue   = [0.45 0.70 1.00];
purple = [0.70 0.50 0.85];

n1 = ceil(npos/2);
n2 = npos - n1 + 1;

cmap_pos = [ ...
    interp1([0 1], [white; blue], linspace(0,1,n1)); ...
    interp1([0 1], [blue; purple], linspace(0,1,n2))];

cmap_pos(n1+1,:) = [];
cmap = [cmap_neg; cmap_pos];

% Tick locations and labels
ticks = 1:3*M;

labels_x = cell(3*M,1);
labels_y = cell(3*M,1);

for k = 1:M
    labels_x{3*k-2} = sprintf('$\\partial z^{%d}$',k);
    labels_x{3*k-1} = sprintf('$\\partial \\beta^{%d}$',k);
    labels_x{3*k}   = sprintf('$\\partial \\alpha^{%d}$',k);

    labels_y{3*k-2} = sprintf('$\\partial F_z^{%d}$',k);
    labels_y{3*k-1} = sprintf('$\\partial T_x^{%d}$',k);
    labels_y{3*k}   = sprintf('$\\partial T_y^{%d}$',k);
end

% Create 2-by-3 tiled layout
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Plot all six panels
ax = gobjects(2,3);

for i = 1:2
    for j = 1:3

        ax(i,j) = nexttile(t,(i-1)*3+j);

        imagesc(ax(i,j),J_all{i,j});
        axis(ax(i,j),'image');

        colormap(ax(i,j),cmap);
        clim(ax(i,j),color_limits);

        xticks(ax(i,j),ticks);

        if i == 2
            xticklabels(ax(i,j),labels_x);
        else
            xticklabels(ax(i,j),{});
        end

        yticks(ax(i,j),ticks);

        if j == 1
            yticklabels(ax(i,j),labels_y);
        else
            yticklabels(ax(i,j),{});
        end

        set(ax(i,j), ...
            'TickLabelInterpreter','latex', ...
            'FontSize',12);

        title(ax(i,j),J_titles{j}, ...
            'Interpreter','latex','FontSize',14);

        hold(ax(i,j),'on');

        for k = 3:3:3*M-1
            xline(ax(i,j),k+0.5,'k-','LineWidth',1);
            yline(ax(i,j),k+0.5,'k-','LineWidth',1);
        end

        hold(ax(i,j),'off');

    end
end

% Shared colorbar along the top
cb = colorbar(ax(1,3),'northoutside');
cb.Layout.Tile = 'north';

% Overall title
title(t,title_string,'Interpreter','latex','FontSize',16);

end