%==========================================================================
clear all
format long
tic
% Build Hamiltonian
basis_num= 10; % number of basis = 2*basis_num + 1
basis = - basis_num:1:basis_num;
kxBZ = - 0.5:0.01:0.5; % 1st BZ
dx = 0.01; xvec = -pi:dx:pi;
band_num = 3;
iter_eps = 1E-6;
%  Physical Parameters
N = 1; % filling number
Delta_c = 12;
U0 = 32.35/N;
kT = 0.01;  Beta = 1.0/kT; % Temperature
vx_nk = zeros(length(xvec),length(kxBZ),band_num);
%%%%%%%%%%%%%%%%%%%%%%%%
Eta = 1.2;  
% ***********************  MAIN  LOOP  ************************************
iter_count = 0;
fprintf('kT=%g, N=%g, Delta_c=%g, U0=%g\n',kT,N,Delta_c,U0);
disp('***************** Start Iteration ! ******************');
alpha = -0.000253425403851 + 0.249876564041172i;
[Eband,Vband] = GetEigens(basis,kxBZ,band_num,U0,alpha,Eta);
figure
subplot(1,2,1)
 plot(kxBZ,Eband(:,1:2),'-o')
 hold on
 idmu = 1;
 %Eta=1.5: 0.920948866662639,  1.095980261883295
 % Eta=1.2, Mu=1.095980263958693;alpha=-0.000253425403851 + 0.249876564041172i
 muvec = 0.80:0.015:1.2;
 Nmu=[];
for Mu = muvec
   [rhox,~] = GetRho(Eband,Vband,xvec,basis,kxBZ,Mu,Beta); %Get the rho(x) in the coor space
   Nmu(idmu) = sum(rhox)*dx;
   % check sum rhox=filling number, to extract new Mu !
   if (mod(iter_count,1)==0)
    fprintf('|alph|=%7.6f,phi=%7.6f*pi,mu=%7.6f,Nmu=%5.3f\n',...
        abs(alpha),phase(alpha)/pi,Mu,Nmu(idmu));
   end
 subplot(1,2,1)
 plot(kxBZ,Mu*ones(length(kxBZ),1),'-')
 subplot(1,2,2)
 plot(muvec(idmu),Nmu(idmu),'-o')
 hold on
   idmu = idmu+1;
end
subplot(1,2,1)
xlabel('k_x (BZ)');ylabel('Band Energy');
subplot(1,2,2)
xlabel('\mu');ylabel('N(\mu)');
hold on




