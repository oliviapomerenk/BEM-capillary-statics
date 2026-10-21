



%% v2

% used ChatGPT

function f = getSurfaceFunctionGUI(L)
% GUI for defining a surface z = f(x,y,xCoG,yCoG,L).
%
% Standard option:
%   f = k*L*((x-xCoG)/L)^2 - k*L*((y-yCoG)/L)^2
%
% Custom option:
%   User-defined anonymous function entered as a string.


    f = [];

    % Create GUI
    fig = uifigure( ...
        'Name', 'Define Surface', ...
        'Position', [200 300 650 470], ...
        'Resize', 'off');

    % Surface selection
    bg = uibuttongroup(fig, ...
        'Title', 'Choose a surface definition', ...
        'Position', [20 75 610 350], ...
        'SelectionChangedFcn', @updateMode);

    % Option 1: Standard surface
        rbStandard = uiradiobutton(bg, ...
            'Text', 'Default surface (saddle)', ...
            'Position', [15 295 560 25]);

        uilabel(bg, ...
            'Position', [35 250 550 35], ...
            'Text', '$$z=kL\left[\left(\frac{x-x_{\mathrm{CoG}}}{L}\right)^2-\left(\frac{y-y_{\mathrm{CoG}}}{L}\right)^2\right]$$', ...
            'Interpreter', 'latex', ...
            'FontSize', 11);

        uilabel(bg, ...
            'Position', [40 210 320 25], ...
            'Text', 'Deformation parameter k (recommend 0 <= k <= 0.1):');

        kField = uieditfield(bg, 'numeric', ...
            'Position', [340 210 70 25], ...
            'Value', 0);

    % Option 2: Custom surface
    rbCustom = uiradiobutton(bg, ...
        'Text', 'Custom user-defined surface', ...
        'Position', [15 165 560 25]);

    uilabel(bg, ...
        'Position', [30 125 560 30], ...
        'Text', 'Enter an expression using x, y, xCoG, yCoG, and L. Please use element-wise operators.');

    funcField = uitextarea(bg, ...
        'Position', [30 30 550 90], ...
        'Value', {''}, ...
        'Enable', 'off');

    % Default selection
    bg.SelectedObject = rbStandard;

    % Buttons
    uibutton(fig, 'push', ...
        'Text', 'OK', ...
        'Position', [390 25 100 30], ...
        'ButtonPushedFcn', @acceptCallback);

    % Wait for user input
    uiwait(fig);

    if isvalid(fig)
        delete(fig);
    end

    % Enable only the selected option's inputs
    function updateMode(~, event)

        if event.NewValue == rbStandard
            kField.Enable = 'on';
            funcField.Enable = 'off';
        else
            kField.Enable = 'off';
            funcField.Enable = 'on';
        end

    end

    % Return the selected function
    function acceptCallback(~, ~)

        if bg.SelectedObject == rbStandard

            % Capture k inside the anonymous function
            kStr = mat2str(kField.Value, 17);
            LStr = mat2str(L, 17);

            f = str2func([ ...
                '@(x,y,xCoG,yCoG) (' kStr ').*(' LStr ').*' ...
                '(((x-xCoG)./(' LStr ')).^2 - ((y-yCoG)./(' LStr ')).^2)' ...
                ]);

        else

            % Read custom expression
            expression = strtrim( ...
                strjoin(string(funcField.Value), ' '));

            if strlength(expression) == 0
                uialert(fig, ...
                    'Please enter a custom surface expression.', ...
                    'Missing Expression');
                return;
            end

            % Construct anonymous function
            try

                % Allowed variables
                allowedVariables = {'x', 'y', 'xCoG', 'yCoG', 'L'};
                
                % Allowed mathematical functions
                allowedFunctions = { ...
                    'sin', 'cos', 'tan', ...
                    'asin', 'acos', 'atan', 'atan2', ...
                    'sinh', 'cosh', 'tanh', ...
                    'asinh', 'acosh', 'atanh', ...
                    'exp', 'log', 'log10', 'sqrt', ...
                    'abs', 'sign', 'floor', 'ceil', 'round', 'fix', ...
                    'mod', 'rem', 'hypot'};
                
                % Extract identifiers from the expression
                exprCheck = regexprep(char(expression), ...
                    '(?<![A-Za-z0-9_.])(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?', '');
                
                identifiers = unique(regexp(exprCheck, '[A-Za-z_]\w*', 'match'));
                
                % Reject names that are neither permitted variables nor functions
                allowedNames = [allowedVariables, allowedFunctions];
                
                invalidNames = identifiers(~ismember(identifiers, allowedNames));
                
                if ~isempty(invalidNames)
                    uialert(fig, ...
                        sprintf(['Unrecognized variable or function: %s\n\n' ...
                        'Use x, y, xCoG, yCoG, L, numbers, ' ...
                        'and supported mathematical functions.'], ...
                        strjoin(invalidNames, ', ')), ...
                        'Invalid Expression');
                    return;
                end

                % Require element-wise operators
                if ~isempty(regexp(char(expression), '(?<!\.)[*/^]', 'once'))
                    uialert(fig, ...
                        ['Use element-wise operators in your expression:' newline ...
                         '  .* for multiplication' newline ...
                         '  ./ for division' newline ...
                         '  .^ for powers'], ...
                        'Invalid Operators');
                    return;
                end

                LStr = mat2str(L, 17);

                expression = regexprep(char(expression), ...
                    '\bL\b', ['(' LStr ')']);
                
                f = str2func(['@(x,y,xCoG,yCoG) ' expression]);

            catch ME
                uialert(fig, ...
                    sprintf('Invalid function expression:\n\n%s', ...
                    ME.message), ...
                    'Invalid Expression');
                return;
            end

        end

        uiresume(fig);
        delete(fig);

    end


end