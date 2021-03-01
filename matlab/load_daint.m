function [ RESDIR ] = load_daint( outfilename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    hostfolder  = outfilename(1:end-8);
    hostfile    = [hostfolder,'/out*'];
    RESDIR      = ['../',outfilename(33:end-8),'/'];
    localfolder = [RESDIR,'.'];
    elafolder   = ['HeLaZ/',outfilename(33:end-8)];
    % Create directory to move the results
    CMD = ['ssh ahoffman@ela.cscs.ch ',...
        'mkdir -p ',elafolder];
    disp(CMD);
    system(CMD);
    
    % Move results file from daint to ela3
    CMD = ['ssh ahoffman@ela.cscs.ch ssh ahoffman@daint.cscs.ch ',...
        'mv ', hostfolder,'/* ',elafolder];
    disp(CMD);
    system(CMD);
    
    % Download results from ela3
    CMD = ['scp -r ahoffman@ela.cscs.ch:',hostfile,' ',localfolder];
    disp(CMD);
    system(CMD);
end

