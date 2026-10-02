function [xcog_stacked, ycog_stacked, L, xmin, xmax, ymin, ymax] = ...
    particlePositionGraphicalGUI(gridSpacing, snapToGrid, M)

    % Developed using ChatGPT

    % Initialize output vectors
    xcog_stacked = [];
    ycog_stacked = [];

    % =========================================================
    % Create figure
    % =========================================================

    fig = uifigure('Name', 'Particle Positions', ...
                   'Position', [100 150 850 650]);

    % =========================================================
    % Input fields
    % =========================================================

    % X limits
    uilabel(fig, ...
        'Position', [30 590 70 22], ...
        'Text', 'X limits:');

    xMinField = uieditfield(fig, 'numeric', ...
        'Position', [100 590 80 22], ...
        'Value', -5);

    xMaxField = uieditfield(fig, 'numeric', ...
        'Position', [190 590 80 22], ...
        'Value', 5);

    % Y limits
    uilabel(fig, ...
        'Position', [30 550 70 22], ...
        'Text', 'Y limits:');

    yMinField = uieditfield(fig, 'numeric', ...
        'Position', [100 550 80 22], ...
        'Value', -5);

    yMaxField = uieditfield(fig, 'numeric', ...
        'Position', [190 550 80 22], ...
        'Value', 5);

    % Characteristic length
    uilabel(fig, ...
        'Position', [30 510 200 22], ...
        'Text', 'Characteristic radius (cm):');

    LField = uieditfield(fig, 'numeric', ...
        'Position', [190 510 80 22], ...
        'Value', 1);

    % Set Grid button
    uibutton(fig, 'push', ...
        'Text', 'Set lims & radius', ...
        'Position', [120 440 100 25], ...
        'ButtonPushedFcn', @setGrid);

    % Particle count display
    countLabel = uilabel(fig, ...
        'Position', [30 390 240 22], ...
        'Text', sprintf('Particles created: 0 / %d', M));

    % Undo button
    uibutton(fig, 'push', ...
        'Text', 'Undo', ...
        'Position', [250 30 100 40], ...
        'ButtonPushedFcn', @undoPoint);

    % Finish button
    uibutton(fig, 'push', ...
        'Text', 'Finish', ...
        'Position', [370 30 100 40], ...
        'ButtonPushedFcn', @finish);

    % =========================================================
    % Axes
    % =========================================================

    ax = uiaxes(fig, ...
        'Position', [320 100 500 450]);

    grid(ax, 'on');
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';

    title(ax, sprintf('Click to place the center of geometry of %d particles.', M));
    xlabel(ax, 'x');
    ylabel(ax, 'y');

    % Initial limits
    xlim(ax, [xMinField.Value xMaxField.Value]);
    ylim(ax, [yMinField.Value yMaxField.Value]);

    % Keep limits fixed
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';

    % Make axes respond to mouse clicks
    ax.ButtonDownFcn = @addPoint;

    % =========================================================
    % Graphics handles and particle counter
    % =========================================================

    particleCenters = gobjects(0);
    particleCircles = gobjects(0);

    % Explicit number of currently existing particles
    nParticles = 0;

    % =========================================================
    % Set Grid
    % =========================================================

    function setGrid(~, ~)

        xmin = xMinField.Value;
        xmax = xMaxField.Value;
        ymin = yMinField.Value;
        ymax = yMaxField.Value;
        L    = LField.Value;

        % Check inputs
        if xmin >= xmax
            uialert(fig, ...
                'The minimum X limit must be smaller than the maximum X limit.', ...
                'Invalid X limits');
            return
        end

        if ymin >= ymax
            uialert(fig, ...
                'The minimum Y limit must be smaller than the maximum Y limit.', ...
                'Invalid Y limits');
            return
        end

        if L <= 0
            uialert(fig, ...
                'Particle radius L must be positive.', ...
                'Invalid particle radius');
            return
        end

        % Set axis limits
        xlim(ax, [xmin xmax]);
        ylim(ax, [ymin ymax]);

        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';

        % Delete existing particle graphics
        if ~isempty(particleCenters)
            validCenters = isvalid(particleCenters);
            delete(particleCenters(validCenters));
        end

        if ~isempty(particleCircles)
            validCircles = isvalid(particleCircles);
            delete(particleCircles(validCircles));
        end

        particleCenters = gobjects(0);
        particleCircles = gobjects(0);

        % Redraw existing particles with new L
        hold(ax, 'on');

        for i = 1:nParticles

            particleCenters(i) = plot(ax, ...
                xcog_stacked(i), ycog_stacked(i), ...
                'ko', ...
                'MarkerFaceColor', 'k', ...
                'MarkerSize', 6);

            particleCircles(i) = rectangle(ax, ...
                'Position', [xcog_stacked(i)-L, ...
                             ycog_stacked(i)-L, ...
                             2*L, ...
                             2*L], ...
                'Curvature', [1 1], ...
                'EdgeColor', 'r', ...
                'LineWidth', 1.5);

        end

        hold(ax, 'off');

        % Update count
        countLabel.Text = sprintf( ...
            'Particles created: %d / %d', nParticles, M);

        % Re-lock axis limits
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';

    end

    % =========================================================
    % Add point
    % =========================================================

    function addPoint(~, ~)

        % Do not allow more than M particles
        if nParticles >= M
            return
        end

        % Save current axis limits
        currentXLim = ax.XLim;
        currentYLim = ax.YLim;

        % Get mouse position
        cp = ax.CurrentPoint;

        if ~snapToGrid
            xClick = cp(1,1);
            yClick = cp(1,2);
        else
            xClick = round(cp(1,1)/gridSpacing) * gridSpacing;
            yClick = round(cp(1,2)/gridSpacing) * gridSpacing;
        end

        % Only accept clicks inside the axes
        if xClick < currentXLim(1) || xClick > currentXLim(2) || ...
           yClick < currentYLim(1) || yClick > currentYLim(2)
            return
        end

        % Store point
        xcog_stacked(end+1,1) = xClick;
        ycog_stacked(end+1,1) = yClick;

        % Increment particle counter
        nParticles = nParticles + 1;

        % Characteristic radius
        L = LField.Value;

        % Keep existing graphics
        hold(ax, 'on');

        % Plot particle center
        particleCenters(nParticles) = plot(ax, ...
            xClick, yClick, ...
            'ko', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 6);

        % Draw particle circle
        particleCircles(nParticles) = rectangle(ax, ...
            'Position', [xClick-L, ...
                         yClick-L, ...
                         2*L, ...
                         2*L], ...
            'Curvature', [1 1], ...
            'EdgeColor', 'r', ...
            'LineWidth', 1.5);

        hold(ax, 'off');

        % Update count
        countLabel.Text = sprintf( ...
            'Particles created: %d / %d', nParticles, M);

        % Restore axis limits
        xlim(ax, currentXLim);
        ylim(ax, currentYLim);

        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';

    end

    % =========================================================
    % Undo
    % =========================================================

    function undoPoint(~, ~)

        % Nothing to undo
        if nParticles == 0
            return
        end

        % Delete graphics for the most recently added particle
        if isvalid(particleCenters(nParticles))
            delete(particleCenters(nParticles));
        end

        if isvalid(particleCircles(nParticles))
            delete(particleCircles(nParticles));
        end

        % Remove graphics handles
        particleCenters(nParticles) = [];
        particleCircles(nParticles) = [];

        % Remove coordinates
        xcog_stacked(nParticles) = [];
        ycog_stacked(nParticles) = [];

        % Decrement particle counter
        nParticles = nParticles - 1;

        % Update count
        countLabel.Text = sprintf( ...
            'Particles created: %d / %d', nParticles, M);

    end

    % =========================================================
    % Finish
    % =========================================================

    function finish(~, ~)

        nPoints = nParticles;

        if nPoints ~= M

            uialert(fig, ...
                sprintf(['You must create exactly %d particles.\n\n' ...
                         'You have currently created %d.'], ...
                         M, nPoints), ...
                'Incorrect Number of Particles');

            return
        end

        % Save final input values
        xmin = xMinField.Value;
        xmax = xMaxField.Value;
        ymin = yMinField.Value;
        ymax = yMaxField.Value;
        L    = LField.Value;

        uiresume(fig);

    end

    % =========================================================
    % Wait for user to finish
    % =========================================================

    uiwait(fig);

    % Delete GUI
    if isvalid(fig)
        delete(fig);
    end

end