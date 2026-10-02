function [x, y, xcm, ycm, zcm, user_set_cm, N, perimeter, Rmin, ds] = drawContourGUI(L, snapToGrid, gridSpacing, M, xcog_stacked, ycog_stacked, dev, Nmin, Nmax, C_BEM, C_curv, thickness)

    % Developed with assistance from ChatGPT

    devFactor = 1 + 5*dev; % wavier particles (in 3D) get more points

    % get approximate minimum gap among particles -- for computing # boundary points
    if M ~= 1

        D = 2*L; % effective particle diameter

        dx = xcog_stacked - xcog_stacked.'; % pairwise center-to-center distances
        dy = ycog_stacked - ycog_stacked.';

        dist_centers = sqrt(dx.^2 + dy.^2);

        dist_centers(dist_centers == 0) = Inf; % ignore self-distances

        min_edge_distance = min(dist_centers(:)) - D; % Minimum edge-to-edge distance

        if min_edge_distance < 0

            disp('Warning: current configuration may involve overlaps.')

            min_edge_distance = input('Assign characteristic inter-particle gap in cm: ');

            disp('Note that an error will be thrown later if particles overlap after the user defines particle shape & gives initial rotations.');

        end

    else

        min_edge_distance = Inf;

    end

    % Check inputs
    if nargin < 1 || ~isscalar(L) || L <= 0
        error('L must be a positive scalar.');
    end

    % Optional inputs
    if nargin < 2
        snapToGrid = false;
    end

    if nargin < 3
        gridSpacing = 0.1;
    end

    % Initialize particle center
    xcm = 0;
    ycm = 0;
    zcm = 0;

    % Flag indicating whether the user explicitly set a displacement
    user_set_cm = 0;

    % =========================================================
    % FIRST STAGE: Choose shape
    % =========================================================

    fig = uifigure('Name', 'Choose Shape', ...
                   'Position', [300 300 450 300]);

    uilabel(fig, ...
        'Position', [50 235 350 25], ...
        'Text', 'Choose the particle shape:', ...
        'FontSize', 14);

    shapeDropdown = uidropdown(fig, ...
        'Items', {'Circle', 'Ellipse', 'Custom'}, ...
        'Value', 'Circle', ...
        'Position', [50 195 200 30], ...
        'ValueChangedFcn', @shapeChanged);

    % Eccentricity input
    eLabel = uilabel(fig, ...
        'Position', [50 145 180 25], ...
        'Text', 'Eccentricity 0 < e < 1:');

    eField = uieditfield(fig, 'numeric', ...
        'Position', [180 145 100 25], ...
        'Value', 0.8);

    % Hide ellipse-specific input initially
    eLabel.Visible = 'off';
    eField.Visible = 'off';

    % Continue button
    uibutton(fig, 'push', ...
        'Text', 'Continue', ...
        'Position', [150 50 150 40], ...
        'ButtonPushedFcn', @chooseShape);

    % Wait for user
    uiwait(fig);

    % Save selected shape
    selectedShape = shapeDropdown.Value;

    % Save eccentricity for ellipse
    if strcmp(selectedShape, 'Ellipse')
        e = eField.Value;
    end

    % Close shape-selection GUI
    if isvalid(fig)
        delete(fig);
    end

    % =========================================================
    % SECOND STAGE: Generate shape or draw custom contour
    % =========================================================

    if strcmp(selectedShape, 'Circle')

        % -----------------------------------------------------
        % Circle of radius L
        % -----------------------------------------------------

        perimeter = 2*pi*L;
        Rmin = L;

        if isfinite(min_edge_distance)
            smallestLengthScale = min(L,min_edge_distance);
        else
            smallestLengthScale = L;
        end

        ds = min(C_BEM*smallestLengthScale, C_curv*Rmin);

        N = ceil(devFactor * perimeter/ds);
        N = min(Nmax,max(N,Nmin));

        theta = linspace(0, 2*pi, N+1)';
        theta(end) = [];

        x = L*cos(theta);
        y = L*sin(theta);

        x = x(:);
        y = y(:);

    elseif strcmp(selectedShape, 'Ellipse')

        % -----------------------------------------------------
        % Ellipse with semimajor axis L and eccentricity e
        % -----------------------------------------------------

        if e <= 0 || e >= 1
            error('The eccentricity must satisfy 0 < e < 1.');
        end

        b = L*sqrt(1-e^2);

        perimeter = pi * (3*(L+b) - sqrt((3*L+b)*(L+3*b))); % Ramanujan's First Approximation

        Rmin = b^2/L;

        if isfinite(min_edge_distance)
            smallestLengthScale = min(L,min_edge_distance);
        else
            smallestLengthScale = L;
        end

        ds = min(C_BEM*smallestLengthScale, C_curv*Rmin);

        N = ceil(devFactor * perimeter/ds);
        N = min(Nmax,max(N,Nmin));

        theta = linspace(0, 2*pi, N+1)';
        theta(end) = [];

        % ellipse
        x = L*cos(theta);
        y = b*sin(theta);

        % irregular shape
        % r = L*(1 + 0.12*cos(theta) + 0.08*cos(2*theta) + 0.05*sin(3*theta));

        % triangle-ish, order 3 rotational symmetry
        % r = L*(1+0.18*cos(3*theta) - 0.03*cos(6*theta));
        %
        % x = r.*cos(theta);
        % y = r.*sin(theta);

        % for a true ellipse, the centroid is at 0, but this is in case the
        % user wants to write in their own function / for testing
        P = polyshape(x,y);

        [xcentroid,ycentroid] = centroid(P);

        x = x - xcentroid;
        y = y - ycentroid;

        x = x(:);
        y = y(:);

    elseif strcmp(selectedShape, 'Custom')

        % -----------------------------------------------------
        % Custom contour
        % -----------------------------------------------------

        fig = uifigure('Name', 'Draw Contour', ...
                       'Position', [100 100 700 600]);

        % Instructions
        uilabel(fig, ...
            'Position', [50 550 600 40], ...
            'Text', ['Draw the contour by clicking and dragging. You do not need to ' ...
                     'return exactly to the starting point; the contour will ' ...
                     'automatically close when you release the mouse.'], ...
            'HorizontalAlignment', 'center', ...
            'WordWrap', 'on', ...
            'FontSize', 11);

        ax = uiaxes(fig, ...
                    'Position', [50 100 600 440]);

        axis(ax, 'equal');

        ax.XGrid = 'on';
        ax.YGrid = 'on';

        title(ax, 'Draw a closed contour');

        % Finish button
        uibutton(fig, 'push', ...
                 'Text', 'Finish', ...
                 'Position', [300 30 100 40], ...
                 'ButtonPushedFcn', @finishDrawing);

        % Draw contour
        h = drawfreehand(ax, ...
                         'Closed', true, ...
                         'Color', 'k');

        h.Waypoints = false(size(h.Position, 1), 1);
        h.InteractionsAllowed = 'none';

        % Wait until user presses Finish
        uiwait(fig);

        % Extract contour coordinates
        pos = h.Position;

        [x,y,N,perimeter,Rmin,ds] = processDrawnParticleContour(pos,L,min_edge_distance, C_BEM, C_curv, Nmin, Nmax, devFactor);

        % Close drawing figure
        if isvalid(fig)
            delete(fig);
        end

    end

    % =========================================================
    % THIRD STAGE: Display normalized contour and select
    % optionally displaced center of mass
    % =========================================================

    fig = uifigure('Name', 'Particle Shape', ...
                   'Position', [100 150 700 600]);

    ax = uiaxes(fig, ...
                'Position', [50 115 600 410]);

    hold(ax, 'on');

    axis(ax, 'equal');

    grid(ax, 'on');

    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';

    % Plot contour
    plot(ax, [x; x(1)], [y; y(1)], ...
         'k-', ...
         'LineWidth', 1.5);

    % Plot nodes
    plot(ax, x, y, ...
         'ko', ...
         'MarkerFaceColor', 'k', ...
         'MarkerSize', 5);

    % Plot center of geometry
    plot(ax, 0, 0, ...
         'r+', ...
         'MarkerSize', 12, ...
         'LineWidth', 2,...
         'HitTest', 'off', ...
         'PickableParts', 'none');

    xlabel(ax, 'x');
    ylabel(ax, 'y');

    title(ax, 'Select a center of mass');

    % ---------------------------------------------------------
    % Homogeneous particle checkbox
    % ---------------------------------------------------------

    homogeneousCheckbox = uicheckbox(fig, ...
        'Position', [50 550 600 25], ...
        'Text', 'Homogeneous particle (CoM computed by geometry)?', ...
        'Value', true, ...
        'FontSize', 20, ...
        'ValueChangedFcn', @homogeneousChanged);

    % Set sensible axis limits
    padding = 0.15 * max([max(abs(x)), max(abs(y)), L]);

    xlim(ax, [min(x)-padding, max(x)+padding]);
    ylim(ax, [min(y)-padding, max(y)+padding]);

    axis(ax, 'equal');

    % ---------------------------------------------------------
    % Instructions for selecting center
    % ---------------------------------------------------------

    uilabel(fig, ...
        'Position', [50 80 600 25], ...
        'Text', ['Click on the plot OR enter numerical displacements. ' ...
                 'Click "Update CoM" to display.'], ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 11);

    % ---------------------------------------------------------
    % x displacement
    % ---------------------------------------------------------

    uilabel(fig, ...
        'Position', [50 48 60 22], ...
        'Text', 'x CoM:');

    xcmField = uieditfield(fig, 'numeric', ...
        'Position', [100 48 80 25], ...
        'Value', 0, ...
        'Enable', 'off');

    % ---------------------------------------------------------
    % y displacement
    % ---------------------------------------------------------

    uilabel(fig, ...
        'Position', [230 48 60 22], ...
        'Text', 'y CoM:');

    ycmField = uieditfield(fig, 'numeric', ...
        'Position', [280 48 80 25], ...
        'Value', 0, ...
        'Enable', 'off');

    % ---------------------------------------------------------
    % z displacement
    % ---------------------------------------------------------

    uilabel(fig, ...
        'Position', [50 12 410 22], ...
        'Text', 'z CoM relative to bottom surface at (xcm,ycm). Default is thickness/2.');

    zcmField = uieditfield(fig, 'numeric', ...
        'Position', [450 12 70 25], ...
        'Value', thickness/2, ...
        'Enable', 'off');

    % ---------------------------------------------------------
    % Update Center button
    % ---------------------------------------------------------

    uibutton(fig, 'push', ...
        'Text', 'Update CoM', ...
        'Position', [420 40 110 40], ...
        'ButtonPushedFcn', @updateCenter);

    % ---------------------------------------------------------
    % Finish button
    % ---------------------------------------------------------

    uibutton(fig, 'push', ...
        'Text', 'Finish', ...
        'Position', [545 40 80 40], ...
        'ButtonPushedFcn', @finishCenter);

    % Allow user to click a point on the axes
    ax.ButtonDownFcn = @selectCenter;

    % Wait for user
    uiwait(fig);

    % Close figure
    if isvalid(fig)
        delete(fig);
    end

    % =========================================================
    % Nested functions
    % =========================================================

    function shapeChanged(~, ~)

        if strcmp(shapeDropdown.Value, 'Ellipse')

            eLabel.Visible = 'on';
            eField.Visible = 'on';

        else

            eLabel.Visible = 'off';
            eField.Visible = 'off';

        end

    end

    function chooseShape(~, ~)

        % Validate ellipse input
        if strcmp(shapeDropdown.Value, 'Ellipse')

            if eField.Value <= 0 || eField.Value >= 1

                uialert(fig, ...
                    'The eccentricity must satisfy 0 < e < 1.', ...
                    'Invalid e');

                return

            end

        end

        uiresume(fig);

    end

    function finishDrawing(~, ~)

        uiresume(fig);

    end

    function homogeneousChanged(~, ~)

        if homogeneousCheckbox.Value

            % Homogeneous particle: CoM is determined by geometry
            xcmField.Enable = 'off';
            ycmField.Enable = 'off';
            zcmField.Enable = 'off';

            user_set_cm = 0;

            xcm = 0;
            ycm = 0;
            zcm = 0;

            xcmField.Value = 0;
            ycmField.Value = 0;

            % Return blue dot to center of geometry
            updateCenterPlot(0, 0);

        else

            % User specifies the CoM
            xcmField.Enable = 'on';
            ycmField.Enable = 'on';
            zcmField.Enable = 'on';

            user_set_cm = 1;

        end

    end

    function selectCenter(~, ~)

        % Do nothing if the particle is homogeneous
        if homogeneousCheckbox.Value
            return
        end

        % Get clicked location
        cp = ax.CurrentPoint;

        if ~snapToGrid

            xcm = cp(1,1);
            ycm = cp(1,2);

        else

            xcm = round(cp(1,1)/gridSpacing)*gridSpacing;
            ycm = round(cp(1,2)/gridSpacing)*gridSpacing;

        end

        % Populate fields
        xcmField.Value = xcm;
        ycmField.Value = ycm;

        % Move blue dot
        updateCenterPlot(xcm, ycm);

    end

    function updateCenter(~, ~)

        % Get values from input fields
        xcm = xcmField.Value;
        ycm = ycmField.Value;

        % Move blue dot
        updateCenterPlot(xcm, ycm);

    end

    function updateCenterPlot(xcm, ycm)

        % Remove previous blue dot
        delete(findobj(ax, 'Tag', 'displacedCenter'));

        % Plot new blue dot
        plot(ax, xcm, ycm, ...
             'bo', ...
             'MarkerFaceColor', 'b', ...
             'MarkerSize', 7, ...
             'Tag', 'displacedCenter');

        drawnow;

    end

    function finishCenter(~, ~)

        % Read final values
        if homogeneousCheckbox.Value

            % Homogeneous particle: CoM is computed from geometry
            xcm = 0;
            ycm = 0;
            zcm = 0;
            user_set_cm = 0;

        else

            % User-specified center of mass
            xcm = xcmField.Value;
            ycm = ycmField.Value;
            zcm = zcmField.Value;
            user_set_cm = 1;

        end

        % Accept current center
        uiresume(fig);

    end

end