function smoothed_y = smooth_curve(x, y, dx)
    % Check if the lengths of x and y are the same
    if length(x) ~= length(y)
        error('Input arrays x and y must have the same length.');
    end
    
    % Initialize the smoothed_y array
    smoothed_y = zeros(size(y));
    
    % Perform smoothing
    for i = 1:length(x)
        indices = find(abs(x - x(i)) <= dx / 2);
        smoothed_y(i) = mean(y(indices));
    end
end