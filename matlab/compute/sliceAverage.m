function [fullAverage, sliceAverages, sliceErrors] = sliceAverage(inputArray, n)
    if mod(length(inputArray), n) ~= 0
        % If not divisible, take a subset of the array
        subsetLength = floor(length(inputArray) / n) * n;
        inputArray = inputArray(1:subsetLength);
    end

    % Reshape the array into n slices
    slicedArray = reshape(inputArray, [], n);

    % Calculate the average over the entire array
    fullAverage = mean(inputArray);

    % Calculate the average for each slice
    sliceAverages = mean(slicedArray);

    % Calculate the error between slice averages
    sliceErrors = abs(sliceAverages - fullAverage);
    
    % % Display results
    % disp(['Full Average: ' num2str(fullAverage)]);
    % for i = 1:n
    %     disp(['Slice ' num2str(i) ' Average: ' num2str(sliceAverages(i))]);
    %     disp(['Error for Slice ' num2str(i) ': ' num2str(sliceErrors(i))]);
    % end
end