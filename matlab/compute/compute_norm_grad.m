function [profile] = compute_norm_grad(profile_folder,xq)
%COMPUTE_NORM_GRAD Compute the normalized gradients L_ref/L_chi=gradchi/chi
    n_stencil = 3;

    ne   = load([profile_folder,'/ne.txt']);
    x    = ne(:,1);
    y    = ne(:,2);
    y    = interp1(x,y,xq);
    % nmvm = ceil(numel(x)*smoothing_ratio);
    % y    = movmean(ne(:,2),nmvm);
    % y    = smooth_curve(x,ne(:,2),dx);
    
    profile.Kne = calculate_derivative(xq,y,n_stencil);
    profile.Kne(:,2) = profile.Kne(:,2)./y;
    profile.ne  = [xq y];

    Te   = load([profile_folder,'/Te.txt']);
    x    = Te(:,1);
    y    = Te(:,2);
    y    = interp1(x,y,xq);
    % nmvm = ceil(numel(x)*smoothing_ratio);
    % y    = movmean(Te(:,2),nmvm);
    % y    = smooth_curve(x,Te(:,2),dx);
    profile.KTe = calculate_derivative(xq,y,n_stencil);
    profile.KTe(:,2) = profile.KTe(:,2)./y;
    profile.Te  = [xq y];

    Ti   = load([profile_folder,'/Ti.txt']);
    x    = Ti(:,1);
    y    = Ti(:,2);
    y    = interp1(x,y,xq);
    % nmvm = ceil(numel(x)*smoothing_ratio);
    % y    = movmean(Ti(:,2),nmvm);
    % y    = smooth_curve(x,Ti(:,2),dx);
    profile.KTi = calculate_derivative(xq,y,n_stencil);
    profile.KTi(:,2) = profile.KTi(:,2)./y;
    profile.Ti  = [xq y];

    % profile.nmvm= nmvm;
    profile.n_stencil= n_stencil;
    profile.xq = xq;
    profile.fold=profile_folder;
end