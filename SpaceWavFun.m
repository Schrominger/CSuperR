function [psi] = SpaceWavFun(kx,kvec,y,x)
% =========================================================================
% Input a wvfun in e^{ikx} basis, convert it into x wave function
% =========================================================================
dim = length(kvec); 
xlength = length(x);
x = reshape(x,xlength,1);
psi = zeros(xlength,1); % complex wavefunction of x psi(x)

  for kid = 1:dim
     psi = psi + y(kid)*exp(1i*(kx+kvec(kid))*x);
  end
  psi = psi/sqrt(2*pi);
end