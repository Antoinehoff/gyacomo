function [variables] = get_miller_GENE_py(prof_folder,rho)

filePath = [prof_folder,'/equilibrium.txt'];

command = ['python3 extract_miller_from_eqdsk.py ',...
            filePath,' ',num2str(rho),' > tmp.txt'];
system(command);

% Specify the path to the output file
filePath = 'tmp.txt';

% Read the content of the file
fileID = fopen(filePath, 'r');
fileContent = textscan(fileID, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fileID);
fileContent = fileContent{1};

% Initialize a structure to store variables
variables = struct();

% Loop through each line in the file
for i = 1:length(fileContent)
    line = strtrim(fileContent{i}); % Remove leading/trailing whitespaces
    
    % Skip empty lines
    if isempty(line)
        continue;
    end
    
    % Split the line into variable name and value
    parts = strsplit(line, '=');
    
    % Skip lines that don't have '=' or have more than two parts
    if length(parts) ~= 2
        continue;
    end
    
    % Extract variable name and value
    variableName = strtrim(parts{1});
    variableValue = str2double(strtrim(parts{2}));
    
    % Check if the conversion to double was successful
    if isnan(variableValue)
        % If conversion fails, try as a string
        variableValue = strtrim(parts{2});
    end
    
    % Store the variable in the structure
    variables.(genvarname(variableName)) = variableValue;
end

system('rm tmp.txt');

end