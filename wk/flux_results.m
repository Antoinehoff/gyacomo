t0 = 100; t1 = 140; [~,it0] = min(abs(t0-Ts2D)); [~,it1] = min(abs(t1-Ts2D)); 
range  = it0:it1;
avg    = mean(Flux_ri(range))
stdev  = std(Flux_ri(range))^(.5)
figure
hist(Flux_ri(range),20)
%%
Gamma_ = [0.62, 0.29];
N_     = [ 256,  128];
L_     = [  50,   25];
P_     = [   2,    2];
J_     = [   1,    1];
NU_    = [0.01, 0.01];
etaB_  = [ 0.5,  0.5];


if 0
%% Fig 3 of Ricci Rogers 2006
fig = figure;
for i = 1:numel(Gamma_)
    semilogy(etaB_(i),Gamma_(i),'o'); hold on;
end
   xlabel('$\eta_B$'); ylabel('$\Gamma^\infty_{part}$') 
end