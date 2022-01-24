function [general]= loadGeneral_GENE(folder,fileNumber,varargin)
% Loads all the input params of a simulation and returns them in a
% structured way, If [show] is provided then results are shown on screen 
% in a compact version.
% To see the way the structure is filled, read the code or try it once.
%
%  INPUTS:
%         folder     -> folder cointaining the simulation data
%         fileNumber -> actual run according to our numeration
%         [show]     -> boolean to disply information on screen
%         [use_h5]   -> define whether h5 files are to be used (default 1) 
%
%       [general]= LOADGENERAL(folder,fileNumber,[show])
%

%%
%we here assume that all the runs are the same in terms of parameters.
fileNumber=fileNumber(1);

use_h5=0;


[allParams]=loadParams(folder,fileNumber,use_h5);

[coord]=loadCoord(folder,fileNumber,allParams.info.x_local,use_h5,allParams.specs.s1.name);

[geometry]=loadGeometry(folder,fileNumber,allParams.info.geom,use_h5,allParams.info.x_local,allParams.box.nx);


% if nargin==3 && varargin{1}==1
%     describe_simulation(folder,fileNumber);
% end

general=struct('coord',coord, 'geometry',geometry, 'folder', folder, 'fileNumber',fileNumber, 'allParams',allParams);

if ~general.allParams.info.x_local       
    for i_s=1:general.allParams.specs.ns
        act_spec=general.allParams.specs.(genvarname(['s',num2str(i_s)])).name;
        if use_h5
            file=fileNameH5(folder,['profiles_',act_spec],fileNumber );
            T=h5read(file,'/temp/T');
            omt=h5read(file,'/temp/omt');
            n=h5read(file,'/density/n');
            omn=h5read(file,'/density/omn');
            x_o_a=h5read(file,'/position/x_o_a');
            x_o_rho_ref=h5read(file,'/position/x_o_rho_ref');
            general.allParams.specs.(genvarname(['s',num2str(i_s)])).profiles=struct('T',T,'n',n,'omn',omn,'omt',omt,'x_o_a',x_o_a,'x_o_rho_ref',x_o_rho_ref);
        else 
            file=fileName_std(folder,['profiles_',act_spec],fileNumber );
            fid=fopen(file);
            nel=numel(regexp(fgetl(fid),'\s+','split'));
            fgetl(fid);
            data=textscan(fid,'%f');
            data=reshape(data{:}, [nel-1 numel(data{:})/(nel-1)]);
            fclose(fid);
            general.allParams.specs.(genvarname(['s',num2str(i_s)])).profiles=...
                struct('T',data(3,:)','n',data(4,:)','omt',data(5,:)','omn',data(6,:)','x_o_a',data(1,:)','x_o_rho_ref',data(2,:)');

        end
    end
      
    general.prefactors=compute_prefactors(general.allParams,general.geometry,numel(coord.ky),numel(coord.z));
end


end
