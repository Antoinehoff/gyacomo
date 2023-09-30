% Get the current directory (wk)
curdir  = pwd;
% Get the scan directory
if 0 % delta kappa scan
    scandir = '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/P2_J1_delta_kappa_scan/';
    % scandir = '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/P4_J2/';
    % scandir = '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/P4_J2_delta_kappa_scan/';
    pnam1 = '$\delta$'; attr1 = 'delta'; pref1 = 0.23;
    pnam2 = '$\kappa$'; attr2 = 'kappa'; pref2 = 1.53;
else % shear safety factor scan
    scandir = '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/P2_J1_PT_sfact_shear_scan/';
    % scandir = '/misc/gyacomo23_outputs/paper_3/DTT_rho85_geom_scan/P2_J1_NT_sfact_shear_scan/';
    pnam1 = '$\hat s$'; attr1 = 'shear'; pref1 = 3.63;
    pnam2 = '$q_0$';    attr2 = 'q0';    pref2 =-2.15;
end
% Get a list of all items in the current directory
contents = dir(scandir);

% give ref value and param names
REFVAL= 0;


% Iterate through the contents
Qxavg = []; Qxerr = []; para1 = []; para2 = [];
for i = 1:length(contents)
    % Check if the item is a directory and not '.' or '..'
    if contents(i).isdir && ~strcmp(contents(i).name, '.') && ~strcmp(contents(i).name, '..')
        % Get and display the name of the subdirectory
        sdubdir = [scandir,contents(i).name];
        disp(['Subdirectory: ' contents(i).name]);
        % Get parameters
        param = read_namelist([sdubdir,'/fort_00.90']);
        para1 = [para1 param.GEOMETRY.(attr1)];
        para2 = [para2 param.GEOMETRY.(attr2)];        
        % Now you are in the subdirectory. You can perform operations here.
        [t_all, Pxi_all, Qxi_all, Pxe_all, Qxe_all] = read_flux_out_XX(sdubdir);
        if(t_all(end) > 40)
            [fullAvg,sliceAvg,sliceErr] = sliceAverage(Qxi_all+Qxe_all,3);
            Qxavg = [Qxavg fullAvg];
            Qxerr = [Qxerr mean(sliceErr)];
        else
            Qxavg = [Qxavg nan];
            Qxerr = [Qxerr nan];
        end
    end
end
%% reshaping, sorting and plotting
p1 = unique(para1);
p2 = unique(para2);
N1 = numel(p1);
N2 = numel(p2);

if para1(1) == para1(2)
    sz = [N2 N1];
    TRANSPOSE = 1;
else
    sz = [N1 N2];
    TRANSPOSE = 0;
end

Zavg = reshape(Qxavg,sz);
Zerr = reshape(Qxerr,sz);
XX   = reshape(para1,sz);
YY   = reshape(para2,sz);

if TRANSPOSE
    Zavg = Zavg';
    Zerr = Zerr';
    XX = XX';
    YY = YY';
end

[~,idx1] = sort(XX(:,1));
[~,idx2] = sort(YY(1,:));
Zavg = Zavg(idx1,idx2);
Zerr = Zerr(idx1,idx2);
XX   = XX(idx1,idx2);
YY   = YY(idx1,idx2);

% compute the 
if REFVAL
    [~,iref1] = min(abs(XX(:,1)-pref1));
    [~,iref2] = min(abs(YY(1,:)-pref2));
    xref  = XX(iref1,iref2);
    yref  = YY(iref1,iref2);
    Qxref = Zavg(iref1,iref2);
    Zavg = (Zavg-Qxref)/Qxref * 100;
    Qxname = '$\langle (Q_{tot}-Q_{ref})/Q_{ref} \rangle_t[\%]$';
    Qrefname = ['$Q_{ref}=$',num2str(Qxref)];
else
    Qxname = '$\langle Q_{tot} \rangle_t$';
end

% Figure
figure
subplot(1,2,1)
% contourf(XX(:,1),YY(1,:),Zavg',13); hold on
[xx_,yy_] = meshgrid(XX(:,1),YY(1,:));
imagesc_custom(xx_,yy_,Zavg'); hold on
if REFVAL
    plot(xref,yref,'xk','MarkerSize',14,'DisplayName',Qrefname)
end
xlabel(pnam1); ylabel(pnam2);
title(Qxname)
colormap(bluewhitered); colorbar;
subplot(1,2,2)
for i = 1:N2
    errorbar(XX(:,i),Zavg(:,i),Zerr(:,i),...
        'DisplayName',[pnam2,'=',num2str(para2(1,i))]);
    hold on;
end
if REFVAL
    plot(xref,0,'xk','MarkerSize',14,'DisplayName',Qrefname)
end
grid on
xlabel(pnam1); ylabel(Qxname);
legend('show');
