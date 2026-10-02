function PlotJacobian(M,J,title_string)



figure;
imagesc(abs(J));
axis image;

% Colormap
colorbar;
colormap(flipud(hot));

% Tick locations: one for each variable within each particle
ticks = 1:3*M;

labels_x = cell(3*M,1);
labels_y = cell(3*M,1);
for k = 1:M
    labels_x{3*k-2} = sprintf('$z^{%d}$',k);
    labels_x{3*k-1} = sprintf('$\\beta^{%d}$',k);
    labels_x{3*k}   = sprintf('$\\alpha^{%d}$',k);
    labels_y{3*k-2} = sprintf('$F_z^{%d}$',k);
    labels_y{3*k-1} = sprintf('$T_x^{%d}$',k);
    labels_y{3*k}   = sprintf('$T_y^{%d}$',k);
end

xticks(ticks);
yticks(ticks);
xticklabels(labels_x);
yticklabels(labels_y);

set(gca, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',12);

xlabel('$\mathbf{q}$','Interpreter','latex','FontSize',16);
ylabel('$\mathbf{R}$','Interpreter','latex','FontSize',16);

title([title_string '$|\mathbf{J}|=|d\mathbf{R}/d\mathbf{q}|$'], ...
    'Interpreter','latex','FontSize',16);

hold on;

for k = 3:3:3*M-1
    xline(k + 0.5,'k-','LineWidth',1);
    yline(k + 0.5,'k-','LineWidth',1);
end

%set(gca,'ColorScale','log');
%clim([10e-1 10e4]);

hold off;

end