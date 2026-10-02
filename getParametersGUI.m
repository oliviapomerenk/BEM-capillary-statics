function [m, thickness, M, is_flipped, snapToGrid, gridSpacing, inputMethod] = getParametersGUI()

    % Developed using ChatGPT

    % Create figure
    fig = uifigure('Name', 'Input Parameters', ...
                   'Position', [500 200 400 425], ...
                   'Resize', 'off');


    %% Number of particles
    uilabel(fig, ...
        'Position', [50 370 120 22], ...
        'Text', 'Number of particles:');

    particleField = uieditfield(fig, 'numeric', ...
        'Position', [180 370 150 22], ...
        'Value', 2);

    %% Mass
    uilabel(fig, ...
        'Position', [50 335 120 22], ...
        'Text', 'Mass (g):');

    massField = uieditfield(fig, 'numeric', ...
        'Position', [180 335 150 22], ...
        'Value', 0.5);

    %% Thickness
    uilabel(fig, ...
        'Position', [50 295 120 22], ...
        'Text', 'Thickness (cm):');

    thicknessField = uieditfield(fig, 'numeric', ...
        'Position', [180 295 150 22], ...
        'Value', 0.25);

    %% Contact line pinning
    uilabel(fig, ...
        'Position', [50 250 200 22], ...
        'Text', 'Contact line pinning:');

    pinningGroup = uibuttongroup(fig, ...
        'Position', [50 170 280 80], ...
        'BorderType', 'none');

    bottomButton = uiradiobutton(pinningGroup, ...
        'Position', [10 40 150 22], ...
        'Text', 'Bottom-pinned');

    topButton = uiradiobutton(pinningGroup, ...
        'Position', [10 10 150 22], ...
        'Text', 'Top-pinned');

    %% Input method
    uilabel(fig, ...
        'Position', [50 145 200 22], ...
        'Text', 'Particle position input method:');

    inputMethodGroup = uibuttongroup(fig, ...
        'Position', [50 90 280 55], ...
        'BorderType', 'none', ...
        'SelectionChangedFcn', @inputMethodChanged);

    graphicalButton = uiradiobutton(inputMethodGroup, ...
        'Position', [10 30 120 22], ...
        'Text', 'Graphical');

    writeInButton = uiradiobutton(inputMethodGroup, ...
        'Position', [10 5 120 22], ...
        'Text', 'Write-in');

    %% Snap to grid
    snapCheckbox = uicheckbox(fig, ...
        'Position', [200 115 150 22], ...
        'Text', 'Snap to grid?', ...
        'Value', false, ...
        'Enable', 'on', ...
        'ValueChangedFcn', @snapChanged);

    %% Grid spacing
    gridSpacingLabel = uilabel(fig, ...
        'Position', [200 80 120 22], ...
        'Text', 'Grid spacing (cm):', ...
        'Enable', 'off');

    gridSpacingField = uieditfield(fig, 'numeric', ...
        'Position', [200 55 150 22], ...
        'Value', 0.1, ...
        'Enable', 'off');

    %% Submit button
    uibutton(fig, 'push', ...
        'Position', [300 15 80 30], ...
        'Text', 'OK', ...
        'ButtonPushedFcn', @submitCallback);

    %% Wait for user input
    uiwait(fig);

    %% Input method callback
    function inputMethodChanged(~, event)

        if event.NewValue == graphicalButton
            % Graphical input: allow grid snapping
            snapCheckbox.Enable = 'on';

            % Grid spacing remains disabled unless snap is checked
            if snapCheckbox.Value
                gridSpacingLabel.Enable = 'on';
                gridSpacingField.Enable = 'on';
            end

        else
            % Write-in input: disable grid controls
            snapCheckbox.Value = false;
            snapCheckbox.Enable = 'off';

            gridSpacingLabel.Enable = 'off';
            gridSpacingField.Enable = 'off';
        end
    end

    %% Snap-to-grid callback
    function snapChanged(~, ~)

        if snapCheckbox.Value
            gridSpacingLabel.Enable = 'on';
            gridSpacingField.Enable = 'on';
        else
            gridSpacingLabel.Enable = 'off';
            gridSpacingField.Enable = 'off';
        end
    end

    %% Submit callback
    function submitCallback(~, ~)

        % Get mass and thickness
        m = massField.Value;
        thickness = thicknessField.Value;
        % Get number of particles
        M = particleField.Value;

        % Check number of particles
        if isempty(M) || ~isfinite(M) || M <= 0 || M ~= round(M)
            uialert(fig, ...
                'Please enter a positive integer number of particles.', ...
                'Invalid Input');
            return
        end

        % Check mass
        if isempty(m) || ~isfinite(m) || m < 0
            uialert(fig, ...
                'Please enter a non-negative mass.', ...
                'Invalid Input');
            return
        end

        % Check thickness
        if isempty(thickness) || ~isfinite(thickness) || thickness < 0
            uialert(fig, ...
                'Please enter a non-negative thickness.', ...
                'Invalid Input');
            return
        end

        % Determine contact line pinning
        if bottomButton.Value
            is_flipped = 0;
        elseif topButton.Value
            is_flipped = 1;
        else
            uialert(fig, ...
                'Please select bottom-pinned or top-pinned.', ...
                'Invalid Input');
            return
        end

        % Determine contour input method
        if graphicalButton.Value
            inputMethod = 'graphical';
        elseif writeInButton.Value
            inputMethod = 'write-in';
        else
            uialert(fig, ...
                'Please select Graphical or Write-in.', ...
                'Invalid Input');
            return
        end

        % Determine snap-to-grid setting
        if strcmp(inputMethod, 'graphical')
            snapToGrid = double(snapCheckbox.Value);
        else
            snapToGrid = 0;
        end

        % Get grid spacing
        if snapToGrid
            gridSpacing = gridSpacingField.Value;

            if isempty(gridSpacing) || ~isfinite(gridSpacing) || gridSpacing <= 0
                uialert(fig, ...
                    'Please enter a positive grid spacing.', ...
                    'Invalid Input');
                return
            end
        else
            gridSpacing = [];
        end

        % Close GUI
        uiresume(fig);
        delete(fig);
    end

end
