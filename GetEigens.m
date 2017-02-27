function [Eband,Vband] = GetEigens(basis,kxvec,band_num,U0,Eta,alpha,Delta_c,N)
%==========================================================================
% Get the static spectrum of H = k^2 +V(x), periodic potential V(x)        
% INPUT:
%     basis:
%     kxvec: kx array;    
%     band_num: the # of bands for output!
%     others: parameters of potential
% OUTPUT:
%     Eband: band energy;    Vband: corresponding eigenvectors
%==========================================================================
dim = length(basis);
kdim = length(kxvec);
Eband = zeros(kdim,band_num);
Vband = zeros(dim,kdim,band_num);

idkx = 1;
for kx = kxvec  %  
       Hk = BuildHam(kx,basis,U0,Eta,alpha,Delta_c,N);
       [V, E] = eig(Hk);
       [eigtmp,ind] = sort(diag(E));
        V = V(:,ind);
       Eband(idkx,1:band_num) = eigtmp(1:band_num);
       % for fixed kx the eigenvectors V
       Vband(:,idkx,1:band_num) =  V(:,1:band_num); 
       idkx = idkx + 1;
end

return
end