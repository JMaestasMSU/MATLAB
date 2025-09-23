% Fixed Point Visualizer with Direct Slope Control and f(x) Display
% This allows direct manipulation of g'(x) at the fixed point and shows both f(x) and g(x)

function slope_controlled_fixed_point_visualizer()
    % Create figure and set up UI with dark mode styling
    fig = figure('Position', [100, 100, 1200, 800], 'Name', 'Fixed Point Visualizer: Complete Control', ...
                'Color', [0.2, 0.2, 0.2]);  % Dark background
    
    % Main plotting area with dark theme
    ax = axes('Position', [0.15, 0.35, 0.7, 0.55], ...
              'Color', [0.15, 0.15, 0.15], ...  % Dark axes background
              'XColor', [0.8, 0.8, 0.8], ...    % Light axis lines
              'YColor', [0.8, 0.8, 0.8]);
    
    % Initialize ALL variables at the top level for nested function access
    fixed_pt = 2.0;
    slope = -1.2;
    f_y_intercept = 1.0;
    x0 = 1.5;
    num_iters = 15;
    func_type = 3;
    curvature = 0.5;
    
    % Variables that will be computed in update_plot
    f = [];
    g = [];
    f_str = '';
    g_str = '';
    
    % Create function type selection buttons
    uicontrol('Style', 'text', 'Position', [50, 300, 150, 20], ...
        'String', 'Function Type:', 'FontWeight', 'bold');
    
    btn1 = uicontrol('Style', 'pushbutton', 'Position', [50, 270, 80, 25], ...
        'String', 'Linear', 'Callback', @(~,~) set_function_type(1));
    btn2 = uicontrol('Style', 'pushbutton', 'Position', [140, 270, 80, 25], ...
        'String', 'Quadratic', 'Callback', @(~,~) set_function_type(2));
    btn3 = uicontrol('Style', 'pushbutton', 'Position', [230, 270, 80, 25], ...
        'String', 'Cubic', 'Callback', @(~,~) set_function_type(3));
    
    % Create all sliders
    slope_slider = uicontrol('Style', 'slider', 'Min', -2.0, 'Max', 2.0, ...
        'Value', slope, 'Position', [50, 240, 300, 20], ...
        'Callback', @update_plot, 'SliderStep', [0.025, 0.1]);
    slope_text = uicontrol('Style', 'text', 'Position', [50, 260, 300, 20], ...
        'String', sprintf('Slope g''(%.1f) = %.1f', fixed_pt, slope), 'Tag', 'slope_text');
    
    fp_slider = uicontrol('Style', 'slider', 'Min', 0.5, 'Max', 4.0, ...
        'Value', fixed_pt, 'Position', [50, 200, 300, 20], ...
        'Callback', @update_plot, 'SliderStep', [0.025, 0.1]);
    fp_text = uicontrol('Style', 'text', 'Position', [50, 220, 300, 20], ...
        'String', sprintf('Fixed point location: %.1f', fixed_pt), 'Tag', 'fp_text');
    
    curve_slider = uicontrol('Style', 'slider', 'Min', -1.0, 'Max', 1.0, ...
        'Value', curvature, 'Position', [50, 160, 300, 20], ...
        'Callback', @update_plot, 'SliderStep', [0.05, 0.2]);
    curve_text = uicontrol('Style', 'text', 'Position', [50, 180, 300, 20], ...
        'String', sprintf('Curvature: %.1f', curvature), 'Tag', 'curve_text');
    
    fyint_slider = uicontrol('Style', 'slider', 'Min', -3.0, 'Max', 3.0, ...
        'Value', f_y_intercept, 'Position', [50, 120, 300, 20], ...
        'Callback', @update_plot, 'SliderStep', [0.033, 0.1]);
    fyint_text = uicontrol('Style', 'text', 'Position', [50, 140, 300, 20], ...
        'String', sprintf('f(x) y-intercept: %.1f', f_y_intercept), 'Tag', 'fyint_text');
    
    x0_slider = uicontrol('Style', 'slider', 'Min', 0.1, 'Max', 5.0, ...
        'Value', x0, 'Position', [50, 80, 300, 20], ...
        'Callback', @update_plot, 'SliderStep', [0.02, 0.1]);
    x0_text = uicontrol('Style', 'text', 'Position', [50, 100, 300, 20], ...
        'String', sprintf('Starting point x0: %.1f', x0), 'Tag', 'x0_text');
    
    iter_slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', 30, ...
        'Value', num_iters, 'Position', [50, 40, 300, 20], ...
        'Callback', @update_plot, 'SliderStep', [0.033, 0.1]);
    iter_text = uicontrol('Style', 'text', 'Position', [50, 60, 300, 20], ...
        'String', sprintf('Iterations: %d', num_iters), 'Tag', 'iter_text');
    
    % Display panels
    equation_panel = uicontrol('Style', 'text', 'Position', [400, 60, 280, 180], ...
        'String', '', 'Tag', 'equation_panel', 'FontSize', 11, ...
        'BackgroundColor', [0.9, 0.9, 1], 'HorizontalAlignment', 'left');
    
    analysis_panel = uicontrol('Style', 'text', 'Position', [700, 60, 280, 180], ...
        'String', '', 'Tag', 'analysis_panel', 'FontSize', 11, ...
        'BackgroundColor', [0.95, 0.95, 0.95], 'HorizontalAlignment', 'left');
    
    theory_panel = uicontrol('Style', 'text', 'Position', [1000, 60, 280, 180], ...
        'String', ['CONVERGENCE THEORY:' char(10) char(10) ...
                  '|g''(p)| < 1: Converges' char(10) ...
                  '|g''(p)| > 1: Diverges' char(10) char(10) ...
                  'g''(p) > 0: Monotonic approach' char(10) ...
                  'g''(p) < 0: Oscillating approach' char(10) char(10) ...
                  'Closer to 0: Faster convergence' char(10) ...
                  'Closer to ±1: Slower convergence' char(10) char(10) ...
                  'GREEN: f(x) - root finding' char(10) ...
                  'BLUE: g(x) - fixed point'], ...
        'FontSize', 10, 'BackgroundColor', [0.9, 0.9, 1], 'HorizontalAlignment', 'left');
    
    % Initial plot
    update_plot();
    
    function set_function_type(type)
        func_type = type;
        
        switch func_type
            case 1
                set(findobj('Tag', 'curve_text'), 'String', 'f(x) slope: 0.2');
                set(curve_slider, 'Enable', 'on');
            case 2
                set(findobj('Tag', 'curve_text'), 'String', sprintf('Curvature: %.1f (Quadratic)', curvature));
                set(curve_slider, 'Enable', 'on');
            case 3
                set(findobj('Tag', 'curve_text'), 'String', sprintf('Curvature: %.1f (Cubic)', curvature));
                set(curve_slider, 'Enable', 'on');
        end
        
        update_plot();
    end
    
    function update_plot(~, ~)
        % Get current values and round to 0.1 increments - ADD DEBUG INFO
        slope = round(get(slope_slider, 'Value') * 10) / 10;
        fixed_pt = round(get(fp_slider, 'Value') * 10) / 10;
        f_y_intercept = round(get(fyint_slider, 'Value') * 10) / 10;
        curvature = round(get(curve_slider, 'Value') * 10) / 10;
        x0 = round(get(x0_slider, 'Value') * 10) / 10;
        num_iters = round(get(iter_slider, 'Value'));
        
        % Debug: Print values to command window
        fprintf('Debug: fixed_pt=%.1f, slope=%.1f, func_type=%d\n', fixed_pt, slope, func_type);
        
        % Update slider positions
        set(slope_slider, 'Value', slope);
        set(fp_slider, 'Value', fixed_pt);
        set(fyint_slider, 'Value', f_y_intercept);
        set(curve_slider, 'Value', curvature);
        set(x0_slider, 'Value', x0);
        
        % Update labels
        set(findobj('Tag', 'slope_text'), 'String', sprintf('Slope g''(%.1f) = %.1f', fixed_pt, slope));
        set(findobj('Tag', 'fp_text'), 'String', sprintf('Fixed point location: %.1f', fixed_pt));
        set(findobj('Tag', 'fyint_text'), 'String', sprintf('f(x) y-intercept: %.1f', f_y_intercept));
        set(findobj('Tag', 'x0_text'), 'String', sprintf('Starting point x0: %.1f', x0));
        set(findobj('Tag', 'iter_text'), 'String', sprintf('Iterations: %d', num_iters));
        
        % Calculate actual root based on function type - DEFINE LOCALLY
        actual_root = fixed_pt;  % Default
        
        % Create f(x) and g(x) based on function type
        switch func_type
            case 1 % Linear
                % For linear: f(x) = m*x + b
                % Use curvature slider to control slope of f(x)
                m = curvature;  % Direct control of f(x) slope
                f = @(x) m * x + f_y_intercept;
                % The root of f(x) is where f(x) = 0, so: m*x + f_y_intercept = 0
                % Therefore: x = -f_y_intercept / m
                if abs(m) > 0.01  % Avoid division by zero
                    actual_root = -f_y_intercept / m;
                end
                % g(x) should have its fixed point at the actual root of f(x)
                g = @(x) actual_root + slope * (x - actual_root);
                f_str = sprintf('f(x) = %.2fx + %.1f', m, f_y_intercept);
                g_str = sprintf('g(x) = %.1f + %.1f(x - %.1f)', actual_root, slope, actual_root);
                set(findobj('Tag', 'curve_text'), 'String', sprintf('f(x) slope: %.1f', curvature));
                
            case 2 % Quadratic
                % For quadratic: f(x) = a(x - h)^2 + k
                % We want f(0) = f_y_intercept and we can control the shape
                % Let's make it so the root is at the fixed_pt slider location
                actual_root = fixed_pt;  
                % f(x) = a(x - actual_root)^2 + c, where f(actual_root) = 0
                % So c = 0, and f(0) = a * actual_root^2 = f_y_intercept
                if abs(actual_root) > 0.01
                    a = f_y_intercept / (actual_root^2) + curvature * 0.1;
                else
                    a = curvature * 0.1;
                end
                f = @(x) a * (x - actual_root).^2;
                g = @(x) actual_root + slope * (x - actual_root) + 0.1 * curvature * (x - actual_root).^2;
                f_str = sprintf('f(x) = %.2f(x - %.1f)²', a, actual_root);
                g_str = sprintf('g(x) = %.1f + %.1f(x - %.1f) + %.2f(x - %.1f)²', ...
                    actual_root, slope, actual_root, 0.1*curvature, actual_root);
                set(findobj('Tag', 'curve_text'), 'String', sprintf('Curvature: %.1f (Quadratic)', curvature));
                
            case 3 % Cubic
                % For cubic: f(x) = a(x - actual_root)^3 + b(x - actual_root)
                % We want f(actual_root) = 0 and can control shape with curvature
                actual_root = fixed_pt;
                a = curvature * 0.1;
                % f(0) = -a*actual_root^3 - b*actual_root = f_y_intercept
                if abs(actual_root) > 0.01
                    b = -(f_y_intercept + a * actual_root^3) / actual_root;
                else
                    b = 0.5;
                end
                f = @(x) a * (x - actual_root).^3 + b * (x - actual_root);
                g = @(x) actual_root + slope * (x - actual_root) + 0.05 * curvature * (x - actual_root).^3;
                f_str = sprintf('f(x) = %.2f(x - %.1f)³ + %.2f(x - %.1f)', ...
                    a, actual_root, b, actual_root);
                g_str = sprintf('g(x) = %.1f + %.1f(x - %.1f) + %.2f(x - %.1f)³', ...
                    actual_root, slope, actual_root, 0.05*curvature, actual_root);
                set(findobj('Tag', 'curve_text'), 'String', sprintf('Curvature: %.1f (Cubic)', curvature));
        end
        
        % Update equation display
        equation_str = sprintf('EQUATIONS:\n\nRoot-finding problem:\n%s = 0\n\nFixed-point form:\n%s\n\nActual root: x = %.2f\nFixed point: x = %.2f\ng''(fixed point) = %.1f', ...
                              f_str, g_str, actual_root, actual_root, slope);
        set(findobj('Tag', 'equation_panel'), 'String', equation_str);
        
        % Plotting
        x_range = [max(0, actual_root - 2.5), actual_root + 2.5];
        x_vals = linspace(x_range(1), x_range(2), 300);
        
        cla(ax);
        hold(ax, 'on');
        
        % Plot f(x)
        try
            f_vals = f(x_vals);
            h1 = plot(ax, x_vals, f_vals, 'g-', 'LineWidth', 2);
        catch
            f_vals = zeros(size(x_vals));
            h1 = plot(ax, x_vals, f_vals, 'g-', 'LineWidth', 2);
        end
        
        % Plot reference lines and main functions
        h2 = plot(ax, x_vals, zeros(size(x_vals)), 'k:', 'LineWidth', 1);
        g_vals = g(x_vals);
        h3 = plot(ax, x_vals, g_vals, 'b-', 'LineWidth', 2);
        h4 = plot(ax, x_vals, x_vals, 'k--', 'LineWidth', 1);
        
        % Plot special points
        h5 = plot(ax, actual_root, actual_root, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        h6 = plot(ax, actual_root, 0, 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
        
        % Tangent line at the correct fixed point location
        tangent_x = [actual_root - 0.8, actual_root + 0.8];
        tangent_y = actual_root + slope * (tangent_x - actual_root);
        h7 = plot(ax, tangent_x, tangent_y, 'r--', 'LineWidth', 2);
        
        % Starting point
        h8 = plot(ax, x0, 0, 'mo', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
        
        % Fixed point iteration (without adding to legend)
        x_current = x0;
        converged = false;
        diverged = false;
        
        colors = jet(num_iters);
        for i = 1:num_iters
            x_next = g(x_current);
            
            if abs(x_next) > 20 || abs(x_next - x_current) > 10
                diverged = true;
                break;
            end
            
            if i > 3 && abs(x_next - x_current) < 1e-8
                converged = true;
            end
            
            % Draw cobweb lines without legend entries
            plot(ax, [x_current, x_current], [x_current, x_next], 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot(ax, [x_current, x_next], [x_next, x_next], 'Color', colors(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot(ax, x_next, x_next, 'o', 'Color', colors(i,:), 'MarkerSize', 5, 'HandleVisibility', 'off');
            
            x_current = x_next;
        end
        
        % Create legend with specific handles only
        legend([h1, h2, h3, h4, h5, h6, h7, h8], ...
               {'f(x)', 'f(x) = 0', 'g(x)', 'y = x', 'Fixed Point', 'Root of f(x)', sprintf('Slope = %.1f', slope), 'Start'}, ...
               'Location', 'best', 'NumColumns', 2);
        
        % Finalize plot
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        title(ax, sprintf('Root Finding vs Fixed Point: f(x)=0 and x=g(x)'));
        grid(ax, 'on');
        xlim(ax, x_range);
        ylim(ax, [min([min(f_vals), min(g_vals), x_range(1)]), max([max(f_vals), max(g_vals), x_range(2)])]);
        
        % Analysis
        analysis_str = sprintf('ANALYSIS:\n\n');
        analysis_str = [analysis_str, sprintf('Root of f(x): x = %.2f\n', actual_root)];
        analysis_str = [analysis_str, sprintf('Fixed point of g(x): x = %.2f\n', actual_root)];
        analysis_str = [analysis_str, sprintf('g''(%.2f) = %.1f\n\n', actual_root, slope)];
        
        if abs(slope) < 1
            if abs(slope) < 0.1
                rate = 'Very Fast';
            elseif abs(slope) < 0.5
                rate = 'Fast';
            else
                rate = 'Moderate';
            end
            analysis_str = [analysis_str, sprintf('CONVERGENT (%s)\n', rate)];
        else
            analysis_str = [analysis_str, 'DIVERGENT\n'];
        end
        
        if converged
            analysis_str = [analysis_str, sprintf('Converged: %.4f\n', x_current)];
        elseif diverged
            analysis_str = [analysis_str, 'Diverged!\n'];
        else
            analysis_str = [analysis_str, sprintf('Current: %.3f\n', x_current)];
        end
        
        set(findobj('Tag', 'analysis_panel'), 'String', analysis_str);
        
        hold(ax, 'off');
    end
end

% Run with: slope_controlled_fixed_point_visualizer()