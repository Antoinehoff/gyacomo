function [ RESDIR ] = load_daint( outfilename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    hostfolder  = outfilename(1:end-8);
    hostfile    = [hostfolder,'/out*'];
    RESDIR      = ['../',outfilename(33:end-8),'/'];
    localfolder = [RESDIR,'.'];
    
    % Download results from ela3
    CMD = ['scp -r ahoffman@daint:',hostfile,' ',localfolder];
    disp(CMD);
    system(CMD);
end

