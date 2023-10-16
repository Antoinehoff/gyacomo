function outputs = calculate_derivative(x, y, stencil_size)
    % Check if the lengths of x and y are the same
    if length(x) ~= length(y)
        error('Input arrays x and y must have the same length.');
    end
    
    % Check if stencil_size is an odd number
    if mod(stencil_size, 2) ~= 1
        error('Stencil size must be an odd number for a central difference stencil.');
    end
    
    % Initialize the derivative array
    derivative = zeros(size(y));
    
    % Compute the derivative using central finite difference
    half_stencil = (stencil_size - 1) / 2;
    for i = (1 + half_stencil):(length(x) - half_stencil)
        x_stencil = x(i - half_stencil : i + half_stencil);
        y_stencil = y(i - half_stencil : i + half_stencil);
        
        % Fit a polynomial and compute the derivative
        p = polyfit(x_stencil, y_stencil, 1);
        derivative(i) = p(1);
    end
    
    % Handle boundary points with forward/backward difference
    derivative(1:half_stencil) = (y(half_stencil + 1) - y(1)) / (x(half_stencil + 1) - x(1));
    derivative(end - half_stencil + 1:end) = (y(end) - y(end - half_stencil)) / (x(end) - x(end - half_stencil));
    outputs = derivative;
end
