function [] = tight_layout_generator(m,n)
% Generate a layout with tight subplots
% m = 3; % number lines
% n = 3; % number of columns
figure
t = tiledlayout(m,n,'TileSpacing','Compact','Padding','Compact');
for i = 1:(m*n)
nexttile
end

end
