function [ RESDIR ] = load_marconi( outfilename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    hostfolder  = outfilename(1:end-8);
    hostfile    = [hostfolder,'/outputs*'];
    RESDIR      = ['../',outfilename(46:end-8),'/'];
    localfolder = [RESDIR,'.'];
    system(['scp -r ahoffman@login.marconi.cineca.it:',hostfile,' ',localfolder])
end

