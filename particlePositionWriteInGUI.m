function [xcog_stacked, ycog_stacked, L] = particlePositionWriteInGUI(M)

    % Developed using ChatGPT

    % Create figure
    fig = uifigure('Name', 'Particle Geometry', ...
                   'Position', [500 150 450 500], ...
                   'Resize', 'off');

    %% Characteristic particle radius
    uilabel(fig, ...
        'Position', [40 455 180 22], ...
        'Text', 'Characteristic radius:');

    radiusField = uieditfield(fig, 'numeric', ...
        'Position', [230 455 150 22], ...
        'Value', 1);

    %% Column labels
    uilabel(fig, ...
        'Position', [115 415 150 22], ...
        'Text', 'x center of geometry');

    uilabel(fig, ...
        'Position', [265 415 150 22], ...
        'Text', 'y center of geometry');

    %% Scrollable panel
    particlePanel = uipanel(fig, ...
        'Position', [25 65 405 340], ...
        'Scrollable', 'on');

    %% Particle center fields
    xFields = gobjects(M,1);
    yFields = gobjects(M,1);

    for i = 1:M

        % Put particle 1 at the top of the content
        yPos = 30*(M-i) + 15;

        % Particle number
        uilabel(particlePanel, ...
            'Position', [25 yPos 50 22], ...
            'Text', sprintf('%d:', i));

        % x center
        xFields(i) = uieditfield(particlePanel, 'numeric', ...
            'Position', [75 yPos 120 22], ...
            'Value', 0);

        % y center
        yFields(i) = uieditfield(particlePanel, 'numeric', ...
            'Position', [225 yPos 120 22], ...
            'Value', 0);
    end

    %% OK button
    uibutton(fig, 'push', ...
        'Position', [350 20 80 30], ...
        'Text', 'OK', ...
        'ButtonPushedFcn', @submitCallback);

    %% Wait for input
    uiwait(fig);

    %% Submit callback
    function submitCallback(~, ~)

        % Get characteristic radius
        L = radiusField.Value;

        if isempty(L) || ~isfinite(L) || L <= 0
            uialert(fig, ...
                'Please enter a positive characteristic radius.', ...
                'Invalid Input');
            return
        end

        % Get particle centers
        xcog_stacked = zeros(M,1);
        ycog_stacked = zeros(M,1);

        for j = 1:M

            xcog_stacked(j) = xFields(j).Value;
            ycog_stacked(j) = yFields(j).Value;

            if ~isfinite(xcog_stacked(j)) || ~isfinite(ycog_stacked(j))
                uialert(fig, ...
                    sprintf('Please enter valid x and y coordinates for particle %d.', j), ...
                    'Invalid Input');
                return
            end
        end

        uiresume(fig);
        delete(fig);
    end

end