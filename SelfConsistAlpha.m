function [ alpha ] = SelfConsistAlpha(basis,kxBZ,xvec,band_num,N,kT,U0,Eta,Delta_c)
%[alpha] = SelfConsistAlpha(basis,kxBZ,xvec,band_num,N,kT,U0,Eta,Delta_c)
%%  Only for 1D now %%===========================================
%   Input the cavity and pump parameters  N,Beta,U0, Eta, Delta_c
%  Output the cavity field: <a> = alpha by a self-consistent solution of
%  the coupled equation of atomic and steady state of cavity light field.
%  Ref. PRL188,073602(2017)
%% **************************************************************************
 alpha  = (rand(1)-0.5)+1i*(rand(1)-0.5);
%alpha = a*exp(1i*b);  % random given initial value of alpha
iter_err = 100000;
iter_count = 0;
iter_eps = 1E-6;
fprintf('kT=%f, Delta_c=%f, U_0*N=%f, N=%f, Eta=%f\n',kT,Delta_c,U0,N,Eta);
disp('***************** Start Iteration ! ******************');
while(iter_err > iter_eps)
     iter_count = iter_count +1;
   [Eband,Vband] = GetEigens(basis,kxBZ,band_num,U0,Eta,alpha,Delta_c,N);
   Mu = FindMu(Eband,Vband,xvec,basis,kxBZ,1/kT,N);% Keep N conserves to find Mu
   rhox = GetRho(Eband,Vband,xvec,basis,kxBZ,Mu,1/kT); %Get the rho(x) in the coor space
%    rho_sb = GetRho(Eband(:,1),Vband(:,:,1),xvec,basis,kxBZ,Mu,1/kT);
%    rho_pb = GetRho(Eband(:,2),Vband(:,:,2),xvec,basis,kxBZ,Mu,1/kT);
%    Nsb = sum(rho_sb)*dx;
%    Npb = sum(rho_pb)*dx;
   % check sum rhox=filling number, to extract new Mu !
   old_alpha = alpha;
   alpha = UpdatePara(old_alpha,rhox,xvec,Eta,U0,Delta_c);
   iter_err = abs(abs(alpha) - abs(old_alpha));
   if (mod(iter_count,10)==0)
    fprintf('Err=%8.6f,alpha =%7.6f %.4f*pi,mu=%6.5f\n',iter_err,...
        abs(alpha),phase(alpha)/pi,Mu);
   end
end
 disp('------------- Fnish Iteration ! ---------');
return
end

