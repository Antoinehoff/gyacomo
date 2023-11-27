function [nuGENE, nuGYAC] = compute_coll(Lref,n_e,T_e,T_i)

    tau    = T_i/T_e;
    
    nuGENE = 2.3031E-5*Lref*(n_e)/(T_e)^2*(24.-log(sqrt(n_e*1.0E13)/T_e*0.001));
    
    nuGYAC = 3/8*sqrt(pi/2)*tau.^(3/2)*nuGENE;

end