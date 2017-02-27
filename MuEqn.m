function [ Fmu ] = MuEqn(Eband,Vband,xvec,basis,kxBZ,Mu,Beta,N)
%MUEQN 此处显示有关此函数的摘要
%   此处显示详细说明
dx = abs(xvec(2) - xvec(1));
rhox = GetRho(Eband,Vband,xvec,basis,kxBZ,Mu,Beta);
Fmu = sum(rhox)*dx - N;
return
end

