function [zcm_height_stacked, alpha_stacked, beta_stacked, gamma_stacked] = initialConditionsGUI( ...
    M, x_nodes_stacked, y_nodes_stacked, xcog_stacked, ycog_stacked, xcm_stacked, ycm_stacked)

    %% Check inputs

    if M < 1 || M > 10 || M ~= floor(M)
        error('M must be an integer between 1 and 10.');
    end

    if ~isequal(size(x_nodes_stacked), size(y_nodes_stacked))
        error('x_nodes_stacked and y_nodes_stacked must be the same size.');
    end

    if size(x_nodes_stacked,1) ~= M
        error('x_nodes_stacked and y_nodes_stacked must have M rows.');
    end

    if ~isequal(size(xcog_stacked), [M 1]) || ...
            ~isequal(size(ycog_stacked), [M 1])
        error('xcog_stacked and ycog_stacked must be M-by-1 vectors.');
    end

    %% Figure dimensions

    rowHeight = 50;

    % Large plot on the right
    plotWidth = 650;
    inputWidth = 500;

    figWidth = inputWidth + plotWidth;
    figHeight = max(500, 150 + M*rowHeight);

    fig = uifigure( ...
        'Name', 'Initial Conditions', ...
        'Position', [200 200 figWidth figHeight], ...
        'Resize', 'off');

    %% Initial condition inputs -- LEFT

    % Column positions
    objectX = 30;
    heightX = 100;
    pitchX = 200;
    rollX = 310;
    yawX = 420;

    % Header
    headerY = figHeight - 60;

    uilabel(fig, ...
        'Position', [objectX headerY 60 25], ...
        'Text', 'Particle index', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    uilabel(fig, ...
        'Position', [heightX headerY 80 25], ...
        'Text', 'Height (z, cm)', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    uilabel(fig, ...
        'Position', [pitchX headerY 90 25], ...
        'Text', 'Pitch (x, deg)', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    uilabel(fig, ...
        'Position', [rollX headerY 90 25], ...
        'Text', 'Roll (y, deg)', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    uilabel(fig, ...
        'Position', [yawX headerY 90 25], ...
        'Text', 'Yaw (z, deg)', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    % Input fields
    heightFields = gobjects(M,1);
    pitchFields = gobjects(M,1);
    rollFields = gobjects(M,1);
    yawFields = gobjects(M,1);

    for i = 1:M

        y = figHeight - 50 - i*rowHeight;

        % Object number
        uilabel(fig, ...
            'Position', [objectX y 60 25], ...
            'Text', sprintf('%d', i), ...
            'HorizontalAlignment', 'center');

        % Height
        heightFields(i) = uieditfield(fig, 'numeric');
        heightFields(i).Position = [heightX y 80 25];
        heightFields(i).Value = 0;

        % Pitch
        pitchFields(i) = uieditfield(fig, 'numeric');
        pitchFields(i).Position = [pitchX y 90 25];
        pitchFields(i).Value = 0;

        % Roll
        rollFields(i) = uieditfield(fig, 'numeric');
        rollFields(i).Position = [rollX y 90 25];
        rollFields(i).Value = 0;

        % Yaw
        yawFields(i) = uieditfield(fig, 'numeric');
        yawFields(i).Position = [yawX y 90 25];
        yawFields(i).Value = 0;

    end

    %% OK button

    uilabel(fig, ...
    'Position', [10 60 550 25], ...
    'Text', 'Vertical translation sets the height of CoM z-coordinate. Rotations are then applied about the CoM.', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'center', ...
    'FontSize', 11);

    uibutton(fig, 'push', ...
        'Position', [250 25 100 35], ...
        'Text', 'OK', ...
        'ButtonPushedFcn', @submitCallback);

    %% Particle plot -- RIGHT

    ax = uiaxes(fig, ...
        'Position', [inputWidth + 20 50 plotWidth - 40 figHeight - 80]);

    hold(ax, 'on');

    % Plot all particle boundary nodes
    scatter(ax, ...
        x_nodes_stacked(:), ...
        y_nodes_stacked(:), ...
        25, ...
        'filled');

    % plot centers of mass
    scatter(ax,xcm_stacked,ycm_stacked,80,'filled');

    % Label each particle at its center of geometry
    for i = 1:M

        text(ax, ...
            xcog_stacked(i), ...
            ycog_stacked(i), ...
            sprintf('%d', i), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 14, ...
            'FontWeight', 'bold');

    end

    axis(ax, 'equal');
    grid(ax, 'on');

    xlabel(ax, 'x');
    ylabel(ax, 'y');

    title(ax, 'Particle Assignments');

    hold(ax, 'off');


    %% Wait for user input

    uiwait(fig);

    %% Submit callback

    function submitCallback(~, ~)

        % Initialize outputs
        zcm_height_stacked = zeros(M,1);
        alpha_stacked = zeros(M,1);
        beta_stacked = zeros(M,1);
        gamma_stacked = zeros(M,1);

        % Retrieve values
        for j = 1:M
            zcm_height_stacked(j) = heightFields(j).Value;
            alpha_stacked(j) = pitchFields(j).Value;
            beta_stacked(j) = rollFields(j).Value;
            gamma_stacked(j) = yawFields(j).Value;
        end

        % Check values
        if any(~isfinite(zcm_height_stacked)) || ...
                any(~isfinite(alpha_stacked)) || ...
                any(~isfinite(beta_stacked)) || ...
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