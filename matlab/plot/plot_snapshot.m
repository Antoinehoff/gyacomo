dir = '/home/ahoffman/HeLaZ/results/CBC/NM_F4_kT_4.5_192x64x24x6x4/';
fname = 'check_phi.out';

filename = [dir,fname];
startRow = 2;

formatSpec = '%24f%2s%f%[^\n\r]';

fileID = fopen(filename,'r');

DIMS = fscanf (fileID,'%d %d %d', [3,1]); % indicate the programme in bar chart “header”

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

dataArray{2} = strtrim(dataArray{2});

fclose(fileID);

field_c = dataArray{:, 1} + 1i * dataArray{:, 3};

% Clear temporary variables
clearvars filename startRow formatSpec fileID dataArray ans;

field_c = reshape(field_c,DIMS');

% Plot the snapshot

field_r = ifourier_GENE(mean(field_c,3));

%
% toplot = abs(fftshift(field_c(:,:,1),2));
toplot = real(fftshift(ifourier_GENE(field_c(:,:,1))));
figure
pclr = pcolor(toplot); set(pclr,'EdgeColor','none');
colormap(bluewhitered); shading interp;