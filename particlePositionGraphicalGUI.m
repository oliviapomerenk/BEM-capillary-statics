



function [xcog_stacked, ycog_stacked, xmin, xmax, ymin, ymax] = ...
    particlePositionGraphicalGUI(gridSpacing, snapToGrid, M, x_nodes_prior, y_nodes_prior, xcm_dist, ycm_dist)

    % Developed using ChatGPT

    % Initialize output vectors
    xcog_stacked = [];
    ycog_stacked = [];

    % =========================================================
    % Create figure
    % =========================================================

    fig = uifigure('Name', 'Particle Positions', ...
                   'Position', [100 150 850 650]);

    % Keep the mouse position updated continuously. On some MATLAB releases,
    % a uifigure's CurrentPoint can be inaccurate over a UIAxes unless a
    % WindowButtonMotionFcn is installed.
    lastMousePosition = [NaN NaN];

    % =========================================================
    % Input fields
    % =========================================================

    % X limits
    uilabel(fig, ...
        'Position', [30 590 70 22], ...
        'Text', 'x limits (cm):');

    xMinField = uieditfield(fig, 'numeric', ...
        'Position', [100 590 80 22], ...
        'Value', -5);

    xMaxField = uieditfield(fig, 'numeric', ...
        'Position', [190 590 80 22], ...
        'Value', 5);

    % Y limits
    uilabel(fig, ...
        'Position', [30 550 70 22], ...
        'Text', 'y limits (cm):');

    yMinField = uieditfield(fig, 'numeric', ...
        'Position', [100 550 80 22], ...
        'Value', -5);

    yMaxField = uieditfield(fig, 'numeric', ...
        'Position', [190 550 80 22], ...
        'Value', 5);

    % % Characteristic length
    % uilabel(fig, ...
    %     'Position', [30 510 200 22], ...
    %     'Text', 'Characteristic radius (cm):');
    % 
    % LField = uieditfield(fig, 'numeric', ...
    %     'Position', [190 510 80 22], ...
    %     'Value', 1);

    % Set Grid button
    uibutton(fig, 'push', ...
        'Text', 'Set display axis limits', ...
        'Position', [120 500 150 25], ...
        'ButtonPushedFcn', @setGrid);

    % Particle count display
    countLabel = uilabel(fig, ...
        'Position', [30 420 260 22], ...
        'Text', sprintf('Particles created: 0 / %d', M),'FontSize',15);

    % Editable table of particle center-of-geometry coordinates.

    uilabel(fig, ...
        'Position', [10 360 200 22], ...
        'Text', 'Enter/change positions here:', ...
        'HorizontalAlignment', 'center');

    % Empty cells allow x and y to be entered separately.
    coordinateData = cell(M, 2);
    coordinateTable = uitable(fig, ...
        'Position', [30 90 260 265], ...
        'Data', coordinateData, ...
        'ColumnName', {'x (cm)', 'y (cm)'}, ...
        'ColumnEditable', [true true], ...
        'ColumnWidth', {100, 100}, ...
        'RowName', arrayfun(@num2str, 1:M, 'UniformOutput', false), ...
        'CellEditCallback', @coordinateEdited);



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
    xlabel(ax, 'x (cm)');
    ylabel(ax, 'y (cm)');

    % Initial limits
    xlim(ax, [xMinField.Value xMaxField.Value]);
    ylim(ax, [yMinField.Value yMaxField.Value]);

    % Keep limits fixed
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';

    % Handle mouse clicks at the figure level and convert the actual mouse
    % pixel location to axes data coordinates. This avoids relying on
    % ax.CurrentPoint, which can be stale after deleting/redrawing all plots.
    fig.Units = 'pixels';
    ax.Units = 'pixels';
    fig.WindowButtonMotionFcn = @trackMousePosition;
    fig.WindowButtonDownFcn = @addPoint;

    % =========================================================
    % Graphics handles and particle counter
    % =========================================================

    particleCenters = gobjects(0);
    particleCircles = gobjects(0);
    particleCOMs = gobjects(0);
    particleLabels = gobjects(0);

    % Explicit number of currently existing particles
    nParticles = 0;

    % =========================================================
    % Redraw one particle's graphics from its stored coordinates
    % =========================================================

    function drawParticle(i)
        particleCenters(i) = plot(ax, ...
            xcog_stacked(i), ycog_stacked(i), ...
            'rx', 'MarkerSize', 8, 'LineWidth', 1.5, ...
            'HitTest', 'off', 'PickableParts', 'none');



        particleCircles(i) = plot(ax, ...
            xcog_stacked(i) + x_nodes_prior, ...
            ycog_stacked(i) + y_nodes_prior, ...
            'Color', 'k', 'LineWidth', 1.5, ...
            'HitTest', 'off', 'PickableParts', 'none');

        particleCOMs(i) = plot(ax, ...
            xcog_stacked(i) + xcm_dist, ...
            ycog_stacked(i) + ycm_dist, ...
            'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6, ...
            'HitTest', 'off', 'PickableParts', 'none');

        % Add the particle index just above and to the right of its center.
        % The offset scales with the current axis limits.
        labelOffsetX = 0.015 * diff(ax.XLim);
        labelOffsetY = 0.015 * diff(ax.YLim);
        particleLabels(i) = text(ax, ...
            xcog_stacked(i) + labelOffsetX, ...
            ycog_stacked(i) + labelOffsetY, ...
            sprintf('%d', i), ...
            'Color', [0 0.25 0.8], ...
            'FontWeight', 'bold', ...
            'FontSize', 10, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'HitTest', 'off', 'PickableParts', 'none');
    end

    function updateCoordinateTable()
        data = cell(M, 2);
        for k = 1:nParticles
            data{k, 1} = xcog_stacked(k);
            data{k, 2} = ycog_stacked(k);
        end
        coordinateTable.Data = data;
        countLabel.Text = sprintf('Particles created: %d / %d', nParticles, M);
    end

    % =========================================================
    % Coordinate table edit callback
    % =========================================================

    function coordinateEdited(src, event)
        if isempty(event.Indices)
            return
        end
        row = event.Indices(1);
        col = event.Indices(2);
        data = src.Data;
        % Explicitly apply the current edit in case the table's Data property
        % has not yet reflected event.NewData when this callback runs.
        data{row, col} = event.NewData;

        % Only allow particles to be entered in sequence, so the output
        % coordinate vectors remain compact and ordered.
        if row > nParticles + 1
            uialert(fig, ...
                'Enter coordinates in the next available row.', ...
                'Particle rows must be sequential');
            updateCoordinateTable();
            return
        end

        [xVal, xStatus] = parseCoordinate(data{row, 1});
        [yVal, yStatus] = parseCoordinate(data{row, 2});

        % A row is incomplete until both coordinates are supplied. The table
        % may return edits as character vectors, even for numeric-looking text.
        if strcmp(xStatus, 'blank') || strcmp(yStatus, 'blank')
            return
        end

        if ~strcmp(xStatus, 'valid') || ~strcmp(yStatus, 'valid')
            uialert(fig, 'Coordinates must be finite numeric values.', ...
                'Invalid coordinates');
            updateCoordinateTable();
            return
        end

        if row == nParticles + 1
            nParticles = nParticles + 1;
            xcog_stacked(nParticles, 1) = xVal;
            ycog_stacked(nParticles, 1) = yVal;
        else
            xcog_stacked(row) = xVal;
            ycog_stacked(row) = yVal;
        end

        redrawAllParticles();
    end

    function [number, status] = parseCoordinate(value)
        % uitable cells can contain [] (untouched), numeric scalars, or text.
        % Normalize all supported forms before validating them.
        number = NaN;
        status = 'invalid';

        if isempty(value)
            status = 'blank';
            return
        end

        if isnumeric(value)
            if ~isscalar(value) || isnan(value)
                status = 'blank';
                return
            end
            if isfinite(value)
                number = double(value);
                status = 'valid';
            end
            return
        end

        if ischar(value)
            value = strtrim(value);
            if isempty(value)
                status = 'blank';
                return
            end
            number = str2double(value);
        elseif isstring(value)
            if ~isscalar(value) || ismissing(value)
                status = 'blank';
                return
            end
            value = strtrim(value);
            if strlength(value) == 0
                status = 'blank';
                return
            end
            number = str2double(value);
        else
            return
        end

        if isscalar(number) && isfinite(number)
            status = 'valid';
        else
            number = NaN;
        end
    end

    function redrawAllParticles()
        hold(ax, 'on');
        deleteValidHandles(particleCenters);
        deleteValidHandles(particleCircles);
        deleteValidHandles(particleCOMs);
        deleteValidHandles(particleLabels);
        particleCenters = gobjects(0);
        particleCircles = gobjects(0);
        particleCOMs = gobjects(0);
        particleLabels = gobjects(0);

        for k = 1:nParticles
            drawParticle(k);
        end
        hold(ax, 'off');
        updateCoordinateTable();
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
    end

    function deleteValidHandles(handles)
        if ~isempty(handles)
            valid = isvalid(handles);
            delete(handles(valid));
        end
    end

    % =========================================================
    % Set Grid
    % =========================================================

    function setGrid(~, ~)

        xmin = xMinField.Value;
        xmax = xMaxField.Value;
        ymin = yMinField.Value;
        ymax = yMaxField.Value;
        % L    = LField.Value;

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

        % if L <= 0
        %     uialert(fig, ...
        %         'Particle radius L must be positive.', ...
        %         'Invalid particle radius');
        %     return
        % end

        % Set axis limits
        xlim(ax, [xmin xmax]);
        ylim(ax, [ymin ymax]);

        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';

        % Redraw all existing particles with the new axis limits.
        redrawAllParticles();

        % Update count
        countLabel.Text = sprintf( ...
            'Particles created: %d / %d', nParticles, M);

        % Re-lock axis limits
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';

    end

    % =========================================================
    % Mouse position tracking
    % =========================================================

    function trackMousePosition(src, ~)
        % Reading CurrentPoint during motion avoids stale/inaccurate pointer
        % coordinates reported by some uifigure/UIAxes combinations.
        lastMousePosition = src.CurrentPoint;
    end

    % =========================================================
    % Add point
    % =========================================================

    function addPoint(~, ~)

        % Do not allow more than M particles
        if nParticles >= M
            return
        end

        % Read the physical mouse location in figure pixels, then map it to
        % the drawable axes rectangle, independent of particle graphics.
        mousePos = lastMousePosition;
        % Fall back to the figure value if no motion event has been received.
        if any(~isfinite(mousePos))
            mousePos = fig.CurrentPoint;
        end
        plotPos = ax.InnerPosition;

        % Ignore clicks outside the actual plotting region (e.g. controls,
        % axis labels, tick labels, and title).
        if numel(mousePos) < 2 || numel(plotPos) < 4 || ...
           plotPos(3) <= 0 || plotPos(4) <= 0 || ...
           mousePos(1) < plotPos(1) || mousePos(1) > plotPos(1) + plotPos(3) || ...
           mousePos(2) < plotPos(2) || mousePos(2) > plotPos(2) + plotPos(4)
            return
        end

        currentXLim = ax.XLim;
        currentYLim = ax.YLim;

        % Convert pixel fractions to data coordinates. Account for reversed
        % axis directions even though the GUI normally uses normal directions.
        xFraction = (mousePos(1) - plotPos(1)) / plotPos(3);
        yFraction = (mousePos(2) - plotPos(2)) / plotPos(4);
        if strcmpi(ax.XDir, 'reverse')
            xFraction = 1 - xFraction;
        end
        if strcmpi(ax.YDir, 'reverse')
            yFraction = 1 - yFraction;
        end

        xClick = currentXLim(1) + xFraction * diff(currentXLim);
        yClick = currentYLim(1) + yFraction * diff(currentYLim);

        if snapToGrid
            xClick = round(xClick/gridSpacing) * gridSpacing;
            yClick = round(yClick/gridSpacing) * gridSpacing;
        end

        % Only accept snapped coordinates that remain inside the axes limits.
        if xClick < currentXLim(1) || xClick > currentXLim(2) || ...
           yClick < currentYLim(1) || yClick > currentYLim(2)
            return
        end

        % Store point
        xcog_stacked(end+1,1) = xClick;
        ycog_stacked(end+1,1) = yClick;

        % Increment particle counter
        nParticles = nParticles + 1;

        % % Characteristic radius
        % L = LField.Value;

        % Plot particle center, contour, and center of mass.
        hold(ax, 'on');
        drawParticle(nParticles);
        hold(ax, 'off');

        % Update coordinate table and particle counter.
        updateCoordinateTable();

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
        % L    = LField.Value;

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









