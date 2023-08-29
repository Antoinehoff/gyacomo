% Load the .fig file
figFile = '/home/ahoffman/paper2_fig/31_CBC_linear_convergence.fig'; % Replace with the actual path to your .fig file
% figFile = 'tmp.fig'; % Replace with the actual path to your .fig file
figHandle = openfig(figFile);

%% Extract data from the figure
% Here, we assume that the plotted data is stored in a single line plot
lineHandle = findobj(figHandle, 'Type', 'line');
xData = get(lineHandle, 'XData');
yData = get(lineHandle, 'YData');
% Extract legend labels
legendHandle = findobj(figHandle, 'Type', 'legend');
legendEntries = flip(get(legendHandle, 'String'));
% legendEntries = legendEntries{end:-1:1};
% Close the loaded figure
close(figHandle);
%%
    figure
for i = numel(xData):-1:1
    % plot(xData{i},yData{i},'DisplayName',legendEntries{i}); hold on
    plot(xData{i},yData{i}./xData{i}.^2,'DisplayName',legendEntries{i}); hold on
end
legend show