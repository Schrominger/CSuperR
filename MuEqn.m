function [ Fmu ] = MuEqn(Eband,Vband,xvec,basis,kxBZ,Mu,Beta,N)
%MUEQN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
dx = abs(xvec(2) - xvec(1));
rhox = GetRho(Eband,Vband,xvec,basis,kxBZ,Mu,Beta);
Fmu = sum(rhox)*dx - N;
return
end

