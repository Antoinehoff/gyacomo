% Define the file path
file_path = [resdir,'fort.42']; % Replace 'your_file.txt' with the actual file path

% Open the file for reading
fileID = fopen(file_path, 'r');

% Check if the file was opened successfully
if fileID == -1
    error('Failed to open the file.');
end

% Initialize an empty array to store the data
data__ = [];

% Read the data from the file line by line
while ~feof(fileID)
    % Read a line from the file
    line = fgetl(fileID);
    
    % Convert the line to a numeric value and append it to the data array
    value = str2double(line);
    
    % Check if the conversion was successful
    if ~isnan(value)
        data__ = [data__; value];
    end
end

% Close the file
fclose(fileID);

% Display the stored data
figure;
plot(data__,'o')
