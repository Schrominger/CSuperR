function [mu_out] = FindMu(Eband,Vband,xvec,basis,kxBZ,Beta,N)
%FINDMU  update new chemical potential for self-consitent calculation
%   此处显示详细说明
%syms Mu Fmu rhox
mu0 = 0.5;
fun = @(Mu)MuEqn(Eband,Vband,xvec,basis,kxBZ,Mu,Beta,N);
mu_out = fzero(fun,mu0);
%fprintf(['*** Chemical Potential = ',num2str(mu_out),' ***\n']);

return
end

