function [ RESDIR ] = load_marconi( outfilename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    hostfolder  = outfilename;
    hostfile    = [hostfolder,'out*'];
    MISCDIR     = ['/misc/HeLaZ_outputs/',outfilename(46:end)];
    RESDIR      = ['../',outfilename(46:end)];
    miscfolder =  [MISCDIR,'.'];
    system(['mkdir -p ',miscfolder]);
    disp(['mkdir -p ',miscfolder]);
    resultfolder = [RESDIR,'.'];
%     CMD = ['rsync -r ahoffman@login.marconi.cineca.it:',hostfolder,'/* ',MISCDIR];
%     disp(CMD);
%     system(CMD);
%     CMD = ['rsync -r --exclude ''outputs_*.h5'' ahoffman@login.marconi.cineca.it:',hostfolder,'/* ',MISCDIR];
%     disp(CMD);
%     system(CMD); 
    % SCP the output file from marconi to misc folder of SPCPC
%     CMD = ['scp -r ahoffman@login.marconi.cineca.it:',hostfile,' ',miscfolder];
    CMD = ['rsync ahoffman@login.marconi.cineca.it:',hostfile,' ',miscfolder];
    disp(CMD);
    system(CMD);
    % Load the fort.90 as well in misc folder
    CMD = ['scp -r ahoffman@login.marconi.cineca.it:',hostfolder,'/fort*.90',' ',miscfolder];
    disp(CMD);
    system(CMD);
    % Put it also in the result directory
    CMD = ['scp -r ahoffman@login.marconi.cineca.it:',hostfolder,'/fort*.90',' ',resultfolder];
    disp(CMD);
    system(CMD);
end
