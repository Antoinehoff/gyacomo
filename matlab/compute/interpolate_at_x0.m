function y0 = interpolate_at_x0(x, y, x0)
    % Check if the lengths of x and y are the same
    if length(x) ~= length(y)
        error('Input arrays x and y must have the same length.');
    end
    
    % Perform linear interpolation using interp1
    y0 = interp1(x, y, x0, 'linear', 'extrap');
end
