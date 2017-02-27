function [ rho ] = GetRho(Eband,Vband,xvec,basis,kxBZ,Mu,Beta)
%GETRHO 
%[ rho ] = GetRho(Eband,Vband,xvec,basis,kxBZ,Mu,Beta)
%   Chemical Potential
%syms Mu
kvec = 1*basis; dk = abs(kxBZ(2)-kxBZ(1));
dim = length(basis); xlen = length(xvec);
band_num = length(Eband(1,:));
xvec = reshape(xvec,xlen,1);
vx_nk = zeros(xlen,length(kxBZ),band_num);

count = zeros(xlen,1);
for idn = 1:band_num
   for kxid = 1:length(kxBZ)
  	e_nk = Eband(kxid,idn);
    v_nk = Vband(:,kxid,idn);
    v_nk = v_nk/norm(v_nk);
    nF=1.0/(1.0+exp(Beta*(e_nk-Mu)));
    vx_nk(:,kxid,idn) = SpaceWavFun(kxBZ(kxid),kvec,v_nk,xvec);% psi_{nk}(x)
    count = count + nF*abs(vx_nk(:,kxid,idn)).^2;
  end
end
  rho = count*dk;
end

