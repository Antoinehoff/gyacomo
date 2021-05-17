function [ RESDIR ] = load_marconi( outfilename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    hostfolder  = outfilename(1:end-8);
    hostfile    = [hostfolder,'/out*'];
    MISCDIR     = ['/misc/HeLaZ_outputs/results/',outfilename(46:end-8),'/'];
    RESDIR      = ['../',outfilename(46:end-8),'/'];
    miscfolder =  [MISCDIR,'.'];
    localfolder = [RESDIR,'.'];
    CMD = ['scp -r ahoffman@login.marconi.cineca.it:',hostfile,' ',miscfolder];
    disp(CMD);
    system(CMD);
    CMD = ['scp -r ahoffman@login.marconi.cineca.it:',hostfolder,'/fort.90',' ',miscfolder];
    disp(CMD);
    system(CMD);
end

