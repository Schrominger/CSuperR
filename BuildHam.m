function H = BuildHam(kx,basis,U0,Eta,alpha,Delta_c,N)
% ================================================================
% kx: quasimomentum
% basis_num: # of basis = 2*num_basis+1
% A, V0: parameters for potential
% Vsl=U0*alpha_m^2*cos^2(x)+2*Eta*alpha_m*cos(\Delta\phi)*cos(x)
% ================================================================
dim =length(basis);
%T = zeros(dim,dim); V1 = zeros(dim,dim); V2 = zeros(dim,dim);
tmp = (kx + 1*basis).^2;
T = diag(tmp);
% U0*
V1 =  0.5*eye(dim) + 0.25*diag(ones(dim-2,1),-2)+ 0.25*diag(ones(dim-2,1),+2);
V1 = U0*abs(alpha)^2*V1;

V2 = 0.5*diag(ones(dim-1,1),-1)+ 0.5*diag(ones(dim-1,1),+1);
V2 = Eta*(alpha+conj(alpha))*V2;
% Delta_c|alpha|^2/N
%V3 = -(Delta_c/N)*(alpha(1))^2*eye(dim);
H = T + V1 + V2 ;%+ V3;

return
end