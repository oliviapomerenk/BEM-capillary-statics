
% v6

function [m, thickness, L, M, is_flipped, snapToGrid, gridSpacing, MassOrDensity, ParticleDensity, N, Nmin, Nmax] = getParametersGUI()

    % Developed using ChatGPT

    % Create figure
    fig = uifigure('Name', 'Input Parameters', ...
                   'Position', [500 150 400 600], ...
                   'Resize', 'off');

    %% Number of particles
    uilabel(fig, ...
        'Position', [50 545 120 22], ...
        'Text', 'Number of particles:');

    particleField = uieditfield(fig, 'numeric', ...
        'Position', [240 545 100 22], ...
        'Value', 2);

    %% Mass or density selection
    uilabel(fig, ...
        'Position', [50 510 140 22], ...
        'Text', 'Input quantity:');

    quantityGroup = uibuttongroup(fig, ...
        'Position', [240 505 170 30], ...
        'BorderType', 'none', ...
        'SelectionChangedFcn', @quantityChanged);

    massButton = uiradiobutton(quantityGroup, ...
        'Position', [0 5 70 22], ...
        'Text', 'Mass', ...
        'Value', true);

    densityButton = uiradiobutton(quantityGroup, ...
        'Position', [80 5 85 22], ...
        'Text', 'Density');

    %% Mass
    massLabel = uilabel(fig, ...
        'Position', [50 470 120 22], ...
        'Text', 'Mass (g):');

    massField = uieditfield(fig, 'numeric', ...
        'Position', [240 470 100 22], ...
        'Value', 0.5);

    %% Density
    densityLabel = uilabel(fig, ...
        'Position', [50 470 120 22], ...
        'Text', 'Density (g/cm^3):', ...
        'Visible', 'off');

    densityField = uieditfield(fig, 'numeric', ...
        'Position', [240 470 100 22], ...
        'Value', 1, ...
        'Visible', 'off');

    %% Thickness
    uilabel(fig, ...
        'Position', [50 435 120 22], ...
        'Text', 'Thickness (cm):');

    thicknessField = uieditfield(fig, 'numeric', ...
        'Position', [240 435 100 22], ...
        'Value', 0.25, ...
        'ValueChangedFcn', @thicknessChanged);

    %% Characteristic length
    uilabel(fig, ...
        'Position', [50 400 190 22], ...
        'Text', 'Radius of bounding circle (cm):');

    lengthField = uieditfield(fig, 'numeric', ...
        'Position', [240 400 100 22], ...
        'Value', 1);

    %% Contact line pinning
    uilabel(fig, ...
        'Position', [50 360 200 22], ...
        'Text', 'Contact line pinning:');

    pinningGroup = uibuttongroup(fig, ...
        'Position', [50 280 280 80], ...
        'BorderType', 'none');

    bottomButton = uiradiobutton(pinningGroup, ...
        'Position', [10 50 150 22], ...
        'Text', 'Bottom-pinned');

    topButton = uiradiobutton(pinningGroup, ...
        'Position', [10 20 150 22], ...
        'Text', 'Top-pinned');

    %% Snap to grid
    snapCheckbox = uicheckbox(fig, ...
        'Position', [50 255 120 22], ...
        'Text', 'Snap to grid?', ...
        'Value', false, ...
        'Enable', 'on', ...
        'ValueChangedFcn', @snapChanged);

    %% Grid spacing
    gridSpacingLabel = uilabel(fig, ...
        'Position', [200 255 120 22], ...
        'Text', 'Grid spacing (cm):', ...
        'Enable', 'off');

    gridSpacingField = uieditfield(fig, 'numeric', ...
        'Position', [320 255 60 22], ...
        'Value', 0.1, ...
        'Enable', 'off');

    %% Boundary points selection
    boundaryGroup = uibuttongroup(fig, ...
        'Position', [50 95 230 135], ...
        'BorderType', 'none', ...
        'SelectionChangedFcn', @boundaryChanged);

    adaptiveButton = uiradiobutton(boundaryGroup, ...
        'Position', [0 120 260 22], ...
        'Text', 'Determine # boundary points adaptively', ...
        'Value', true);

    assignButton = uiradiobutton(boundaryGroup, ...
        'Position', [0 15 230 22], ...
        'Text', 'Assign # boundary points:');

    %% Adaptive point-count limits
    minPointsLabel = uilabel(fig, ...
        'Position', [180 185 120 22], ...
        'Text', 'Min # points:');

    minPointsField = uieditfield(fig, 'numeric', ...
        'Position', [280 185 100 22], ...
        'Value', 100, ...
        'Enable', 'on');

    maxPointsLabel = uilabel(fig, ...
        'Position', [180 155 120 22], ...
        'Text', 'Max # points:');

    maxPointsField = uieditfield(fig, 'numeric', ...
        'Position', [280 155 100 22], ...
        'Value', 400, ...
        'Enable', 'on');

    %% Fixed number of boundary points
    boundaryField = uieditfield(fig, 'numeric', ...
        'Position', [280 110 100 22], ...
        'Value', 100, ...
        'Enable', 'off');

    %% Submit button
    uibutton(fig, 'push', ...
        'Position', [300 20 80 30], ...
        'Text', 'OK', ...
        'ButtonPushedFcn', @submitCallback);

    %% Wait for user input
    uiwait(fig);

    %% Quantity selection callback
    function quantityChanged(~, event)

        if event.NewValue == massButton
            massLabel.Visible = 'on';
            massField.Visible = 'on';

            densityLabel.Visible = 'off';
            densityField.Visible = 'off';
        else
            massLabel.Visible = 'off';
            massField.Visible = 'off';

            densityLabel.Visible = 'on';
            densityField.Visible = 'on';
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

    %% Boundary points callback
    function boundaryChanged(~, event)

        if event.NewValue == assignButton
            boundaryField.Enable = 'on';
            minPointsField.Enable = 'off';
            maxPointsField.Enable = 'off';
        else
            boundaryField.Enable = 'off';
            minPointsField.Enable = 'on';
            maxPointsField.Enable = 'on';
        end
    end

    %% Thickness callback
    function thicknessChanged(~, ~)

        if thicknessField.Value == 0
            densityLabel.Text = 'Density (g/cm^2):';
        else
            densityLabel.Text = 'Density (g/cm^3):';
        end
    end

    %% Submit callback
    function submitCallback(~, ~)

        % Get number of particles
        M = particleField.Value;

        if isempty(M) || ~isfinite(M) || M <= 0 || M ~= round(M)
            uialert(fig, ...
                'Please enter a positive integer number of particles.', ...
                'Invalid Input');
            return
        end

        % Get mass or density
        if massButton.Value
            MassOrDensity = 1;
            m = massField.Value;
            ParticleDensity = NaN;

            if isempty(m) || ~isfinite(m) || m < 0
                uialert(fig, ...
                    'Please enter a non-negative mass.', ...
                    'Invalid Input');
                return
            end
        else
            MassOrDensity = 0;
            m = NaN;
            ParticleDensity = densityField.Value;

            if isempty(ParticleDensity) || ...
                    ~isfinite(ParticleDensity) || ParticleDensity < 0
                uialert(fig, ...
                    'Please enter a non-negative density.', ...
                    'Invalid Input');
                return
            end
        end

        % Get thickness and characteristic length
        thickness = thicknessField.Value;
        L = lengthField.Value;

        if isempty(thickness) || ~isfinite(thickness) || thickness < 0
            uialert(fig, ...
                'Please enter a non-negative thickness.', ...
                'Invalid Input');
            return
        end

        if isempty(L) || ~isfinite(L) || L <= 0
            uialert(fig, ...
                'Please enter a positive characteristic length.', ...
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

        snapToGrid = double(snapCheckbox.Value);

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

        % Get boundary point settings
        if assignButton.Value
            N = boundaryField.Value;

            if isempty(N) || ~isfinite(N) || N <= 0 || N ~= round(N)
                uialert(fig, ...
                    'Please enter a positive integer number of boundary points.', ...
                    'Invalid Input');
                return
            end

            % Adaptive limits are unused in fixed-N mode
            Nmin = NaN;
            Nmax = NaN;

        else
            % Adaptive point-count limits
            N = NaN;
            Nmin = minPointsField.Value;
            Nmax = maxPointsField.Value;

            if isempty(Nmin) || ~isfinite(Nmin) || ...
                    Nmin <= 0 || Nmin ~= round(Nmin)
                uialert(fig, ...
                    'Please enter a positive integer minimum number of points.', ...
                    'Invalid Input');
                return
            end

            if isempty(Nmax) || ~isfinite(Nmax) || ...
                    Nmax <= 0 || Nmax ~= round(Nmax)
                uialert(fig, ...
                    'Please enter a positive integer maximum number of points.', ...
                    'Invalid Input');
                return
            end

            if Nmin > Nmax
                uialert(fig, ...
                    'Min # points must be less than or equal to Max # points.', ...
                    'Invalid Input');
                return
            end
        end

        % Close GUI
        uiresume(fig);
        delete(fig);
    end

end