


%% v3



function [zcm_height_stacked, gamma_stacked] = initialConditionsGUI( ...
    M, x_nodes_stacked, y_nodes_stacked, xcog_stacked, ycog_stacked, ...
    xcm_stacked, ycm_stacked, L)

    %% Check inputs

    if ~isequal(size(x_nodes_stacked), size(y_nodes_stacked))
        error('x_nodes_stacked and y_nodes_stacked must be the same size.');
    end

    if size(x_nodes_stacked, 1) ~= M
        error('x_nodes_stacked and y_nodes_stacked must have M rows.');
    end

    if ~isequal(size(xcog_stacked), [M 1]) || ...
       ~isequal(size(ycog_stacked), [M 1])
        error('xcog_stacked and ycog_stacked must be M-by-1 vectors.');
    end

    if ~isequal(size(xcm_stacked), [M 1]) || ...
       ~isequal(size(ycm_stacked), [M 1])
        error('xcm_stacked and ycm_stacked must be M-by-1 vectors.');
    end

    %% Figure dimensions

    rowHeight = 50;

    plotWidth = 650;
    inputWidth = 300;

    figWidth = inputWidth + plotWidth;
    figHeight = 650;

    fig = uifigure( ...
        'Name', 'Initial Conditions', ...
        'Position', [200 200 figWidth figHeight], ...
        'Resize', 'off');

    %% Scrollable input panel -- LEFT
    panelLeft = 10;
    panelBottom = 135;
    panelWidth = inputWidth - 20;
    panelHeight = figHeight - 70 - panelBottom;

    inputPanel = uipanel(fig, ...
        'Position', [panelLeft panelBottom panelWidth panelHeight], ...
        'Scrollable', 'on', ...
        'BorderType', 'none');

    % Total height of the contents inside the panel
    contentHeight = 40 + M*rowHeight;

    % Column positions within the panel
    objectX = 10;
    heightX = 80;
    yawX = 180;

    %% Table headers

    headerY = contentHeight - 30;

    uilabel(inputPanel, ...
        'Position', [objectX headerY 60 25], ...
        'Text', 'Particle', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    uilabel(inputPanel, ...
        'Position', [heightX headerY 90 25], ...
        'Text', 'Height (z, cm)', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    uilabel(inputPanel, ...
        'Position', [yawX headerY 90 25], ...
        'Text', 'Yaw (z, deg)', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    %% Input fields

    heightFields = gobjects(M, 1);
    yawFields = gobjects(M, 1);
    particleLabels = gobjects(M,1);

    for i = 1:M

        y = contentHeight - 30 - i*rowHeight;

        % Particle number
        uilabel(inputPanel, ...
            'Position', [objectX y 60 25], ...
            'Text', sprintf('%d', i), ...
            'HorizontalAlignment', 'center');

        % Height
        heightFields(i) = uieditfield(inputPanel, 'numeric');
        heightFields(i).Position = [heightX y 80 25];
        heightFields(i).Value = 0;

        % Yaw
        yawFields(i) = uieditfield(inputPanel, 'numeric');
        yawFields(i).Position = [yawX y 90 25];
        yawFields(i).Value = 0;
        yawFields(i).ValueChangedFcn = ...
            @(src, event) updateParticlePlot(i);

    end

    %% Instructions and OK button -- FIXED

    uilabel(fig, ...
        'Position', [25 55 inputWidth-50 65], ...
        'Text', ['Vertical translation sets the height of the ' ...
                 'CoM z-coordinate. Yaw rotation is applied ' ...
                 'about the CoM (blue).'], ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'center', ...
        'WordWrap', 'on', ...
        'FontSize', 12);

    uibutton(fig, 'push', ...
        'Position', [120 15 100 35], ...
        'Text', 'OK', ...
        'ButtonPushedFcn', @submitCallback);

    %% Particle plot -- RIGHT

    ax = uiaxes(fig, ...
        'Position', [inputWidth + 20 50 plotWidth - 40 figHeight - 80]);

    hold(ax, 'on');

    % Plot each particle's boundary nodes separately
    for i = 1:M
        scatter(ax, ...
            x_nodes_stacked(i,:), ...
            y_nodes_stacked(i,:), ...
            15, 'k', 'filled', ...
            'Tag', sprintf('particle%d', i));
    end

    % Plot centers of mass
    scatter(ax, xcm_stacked, ycm_stacked, 80, 'blue', 'filled');

    % Label each particle at its center of geometry
    for i = 1:M
        particleLabels(i) = text(ax, ...
            xcog_stacked(i)+L/5, ...
            ycog_stacked(i)+L/5, ...
            sprintf('%d', i), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 14, ...
            'FontWeight', 'bold');
    end

    axis(ax, 'equal');
    grid(ax, 'on');

    xlabel(ax, 'x (cm)');
    ylabel(ax, 'y (cm)');
    title(ax, 'Particle Assignments');

    uilabel(fig, ...
    'Position', [inputWidth + 20 10 plotWidth - 40 30], ...
    'Text', 'If particles overlap, an error will be thrown.', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 12, ...
    'FontColor', [0.8 0 0]);

    hold(ax, 'off');

    %% Wait for user input

    uiwait(fig);

    %% Nested functions

    function updateParticlePlot(i)

        % Current yaw angle in radians
        theta = deg2rad(yawFields(i).Value);

        % Original boundary nodes
        x0 = x_nodes_stacked(i,:);
        y0 = y_nodes_stacked(i,:);

        % Center of mass
        xc = xcm_stacked(i);
        yc = ycm_stacked(i);

        % Counter-clockwise rotation in the x-y plane
        x_rot = xc + (x0 - xc)*cos(theta) - (y0 - yc)*sin(theta);
        y_rot = yc + (x0 - xc)*sin(theta) + (y0 - yc)*cos(theta);

        % Find and update this particle's scatter object
        h = findobj(ax, ...
            'Type', 'Scatter', ...
            'Tag', sprintf('particle%d', i));

        if ~isempty(h)
            h.XData = x_rot;
            h.YData = y_rot;
        end

        % Rotate the particle label's position about the center of mass
        x_label0 = xcog_stacked(i) + L/5;
        y_label0 = ycog_stacked(i) + L/5;
        
        particleLabels(i).Position(1) = ...
            xc + (x_label0 - xc)*cos(theta) - (y_label0 - yc)*sin(theta);
        
        particleLabels(i).Position(2) = ...
            yc + (x_label0 - xc)*sin(theta) + (y_label0 - yc)*cos(theta);

        drawnow;

    end

    function submitCallback(~, ~)

        % Initialize outputs
        zcm_height_stacked = zeros(M, 1);
        gamma_stacked = zeros(M, 1);

        % Retrieve values
        for j = 1:M
            zcm_height_stacked(j) = heightFields(j).Value;
            gamma_stacked(j) = yawFields(j).Value;
        end

        % Check values
        if any(~isfinite(zcm_height_stacked)) || ...
           any(~isfinite(gamma_stacked))

            uialert(fig, ...
                'All values must be finite numbers.', ...
                'Invalid Input');

            return
        end

        % Close GUI
        uiresume(fig);
        delete(fig);

    end

end