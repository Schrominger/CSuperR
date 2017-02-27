function [ new_alpha ] = UpdatePara(~,rho,xvec,Eta,U0,Delta_c)
%UPDATEPARA 
%[ new_alpha ] = UpdatePara(old_alpha,rho,xvec,Eta,U0,Delta_c)
%   finished !
% xvec: spatial coor
% array -- rho in space coordinate
dx =abs(xvec(2)- xvec(1));
new_alpha = (Eta*sum((rho.').*cos(xvec))*dx)/(Delta_c - sum(U0*(rho.').*cos(xvec).^2)*dx);
%% History
% DESCRIPTIVE TEXT
%20170224
  % a1 = sum((rho.').*(Eta*cos(xvec) + U0*real(old_alpha)*cos(xvec).^2))*dx;
  % a1 = a1/Delta_c;
  % 
  % a2 = sum((U0*imag(old_alpha))*(rho.').*(cos(xvec).^2))*dx;
  % a2 = a2/Delta_c;
  %new_alpha = a1+1i*a2
%20170225
  % phi = old_alpha(2);
  % tmp = Eta*cos(phi)*cos(xvec) + U0*(old_alpha(1))*cos(xvec).^2;
  % tmp = (rho.').*tmp;
  % new_alpha(1)= sum(tmp)*dx/Delta_c;
  % new_alpha(2) = old_alpha(2);
%20170226: update modulus of alpha
%     tmp = Eta*cos(phi)*sum((rho.').*cos(xvec))*dx;
%     tmp2 = Delta_c - sum(U0*(rho.').*cos(xvec).^2)*dx;
%     new_alpha(1) = tmp/tmp2;
%     new_alpha(2) = old_alpha(2);

  
return
end

