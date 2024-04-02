function [out] = read_flux_out_XX(folderPath,PLOT,nmvm)
    % Check if the prompt_string is provided as an argument
    if nargin < 1
        % If not provided, prompt the user for input
        prompt_string = 'Enter the path to the simulation dir: ';
        % Get the input from the user
        folderPath = input(prompt_string, 's');   
        PLOT = 1;
        nmvm = 1;
    elseif nargin == 1
            PLOT = 0;
            nmvm = 1;
    elseif nargin == 2
            nmvm = 1;
    end

    % Initialize empty arrays to store the values
    t_all   = [];
    Pxi_all = [];
    Qxi_all = [];
    Pxe_all = [];
    Qxe_all = [];

    % List all files in the folder with the pattern out_XX.txt
    filePattern = fullfile(folderPath, 'out_*.txt');
    files = dir(filePattern);

    % Regular expression pattern to match numerical values
    pattern = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';

    % Loop through each file
    for i = 1:length(files)
        filename = fullfile(folderPath, files(i).name);
        fid = fopen(filename, 'r');

        if fid == -1
            fprintf('Error opening file: %s\n', filename);
            continue;
        end

        % Loop through each line of the file
        while ~feof(fid)
            line = fgetl(fid);

            % Check if the line starts with "|t"
            if startsWith(line, '|t')
                % Check if the line contains numerical values
                matches = regexp(line, pattern, 'match');

                % Convert the matched strings to numerical values
                values = str2double(matches);
                % If matches are found, extract the numerical values and store them in the arrays
                if ~isempty(matches)
                    values = str2double(matches);
                    if numel(values) == 4
                        t_all   = [t_all, values(1)];
                        Pxi_all = [Pxi_all, values(3)];
                        Qxi_all = [Qxi_all, values(4)];
                    elseif numel(values) == 3
                        t_all   = [t_all, values(1)];
                        Pxi_all = [Pxi_all, values(2)];
                        Qxi_all = [Qxi_all, values(3)];
                    elseif numel(values) == 5
                        t_all   = [t_all,   values(1)];
                        Pxi_all = [Pxi_all, values(2)];
                        Qxi_all = [Qxi_all, values(3)];
                        Pxe_all = [Pxe_all, values(4)];
                        Qxe_all = [Qxe_all, values(5)];
                    end
                end
            end
        end

        % Close the file
        fclose(fid);
    end

    if PLOT
    figure
    subplot(211)
    x_ = movmean(t_all,nmvm);
    y_ = movmean(Pxi_all,nmvm);
    plot(x_,y_,'r','DisplayName','ions'); hold on;
    if(numel(t_all)==numel(Pxe_all))
        y_ = movmean(Pxe_all,nmvm);
        plot(x_,y_,'-.b','DisplayName','electrons'); hold on;
    end
    xlabel('$tc_s/R$'); ylabel('$\Gamma_{x}$');
    legend('show')
    title('Radial particle flux')
    subplot(212)
    x_ = movmean(t_all,nmvm);
    y_ = movmean(Qxi_all,nmvm);
    plot(x_,y_,'r','DisplayName','ions'); hold on;
    if(numel(t_all)==numel(Qxe_all))
        y_ = movmean(Qxe_all,nmvm);
        plot(x_,y_,'-.b','DisplayName','electrons'); hold on;
    end
    xlabel('$tc_s/R$'); ylabel('$Q_{x}$');
    legend('show');
    title('Radial heat flux')
    end


    out.t   = t_all;
    out.Pxi = Pxi_all;
    out.Qxi = Qxi_all;
    out.Pxe = Pxe_all;
    out.Qxe = Qxe_all;
end
