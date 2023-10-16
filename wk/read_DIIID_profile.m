function data = read_DIIID_profile(filePath)
    if exist([filePath,'.mat'])
        tmp = load([filePath,'.mat']);
        data = tmp.data;
    else
    % Read the content of the file
    fileID = fopen(filePath, 'r');
    fileContent = textscan(fileID, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fileID);
    fileContent = fileContent{1};
    
    % Initialize variables
    data = struct();
    
    % Loop through each line in the file
    for i = 1:length(fileContent)
        line = strtrim(fileContent{i}); % Remove leading/trailing whitespaces
        
        % Check if the line starts with '201'
        if startsWith(line, '201')
            % Extract column titles with potential units
            titles = strsplit(line, ' ');
            titles(cellfun(@isempty, titles)) = [];  % Remove empty cells
            columnTitle = titles{3};  % Assuming the data starts from the 3rd column
            
            % Remove potential units from the column title
            cleanTitle = regexprep(columnTitle, '\(\S*\)', ''); % Remove everything inside parentheses
            
            % Initialize data field for the current variable
            data.(genvarname(cleanTitle)).x     = [];
            data.(genvarname(cleanTitle)).y     = [];
            data.(genvarname(cleanTitle)).K     = [];
            
            % Extract data for the current column
            data.(genvarname(cleanTitle)).psinorm = sscanf(line, '%f %f %f', [3, inf])';
        else
            % Extract data for the current section
            currentSectionData = sscanf(line, '%f %f %f');
            data.(genvarname(cleanTitle)).x = [data.(genvarname(cleanTitle)).x; currentSectionData(1)];
            data.(genvarname(cleanTitle)).y = [data.(genvarname(cleanTitle)).y; currentSectionData(2)];
            data.(genvarname(cleanTitle)).K = [data.(genvarname(cleanTitle)).K;-currentSectionData(3)];
        end
    end
    % unit conversions
    data.ne.y = data.ne.y * 10; % convert to 10^19 m-3
    save([filePath,'.mat'], 'data');
    end
end
