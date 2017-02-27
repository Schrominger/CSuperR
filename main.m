%% Iteration calculation
%=======================2017.02.22=========================================
%1). initialize alpha=(a1,a2)=a1+i*a2
%    V(x)=U0*abs(alpha)^2*cos^2(x)+2*Eta*abs(alpha)*cos(phase(alpha))*cos(x)
%2). Diagonalize H = k^2 + V(x) => eigenenegies En(k) and psi_{n,k}(x)in BZ
%3). calculate rho(x) by En(k) and psi_{n,k} in BZ
%4). Using rho(x) to calculate the new alpha=a1+i*a2,then RETURN TO step i)
%==========================================================================

clear all
format long
tic
% Build Hamiltonian
basis_num= 10; % number of basis = 2*basis_num + 1
basis = - basis_num:1:basis_num;
kxBZ = - 0.5:0.01:0.5; % 1st BZ
%xvec = linspace(0,2*pi,length(kxBZ));
%dx = xvec(2) - xvec(1);
dx = 0.05; xvec = 0:dx:2*pi;
band_num = 4;
iter_eps = 1E-8;
%  Physical Parameters
N = 1.0; % filling number
Delta_c = 12;
U0 = 32.35/N;
kT = 0.01;  Beta = 1.0/kT; % Temperature
vx_nk = zeros(length(xvec),length(kxBZ),band_num);
%*********  Input Variables   **********
EtaVec = 0.8:0.05:1.4;
% ***********************  MAIN  LOOP  ************************************
alpha_vec = zeros(length(EtaVec),1);
%parpool(2)
parfor idE = 1:length(EtaVec)
alpha_vec(idE) = SelfConsistAlpha(basis,kxBZ,xvec,band_num,N,kT,U0,EtaVec(idE),Delta_c)
end
%%
save('fig2.mat');
figure;
plot(EtaVec,abs(alpha_vec),'k-o')
%% Plot %%%%%%%%%%%%%%%%%%%%
figure;
subplot(1,2,1)
plot(kxBZ,Eband(:,1:2),'-',kxBZ,Mu*ones(length(kxBZ)),'r-');
hold on;
xlabel('k_x (BZ)'); ylabel('Band Energy');

Vsl1 = Eta*(alpha+conj(alpha))*cos(xvec);
Vsl2 = U0*abs(alpha)^2*cos(xvec).^2;
subplot(1,2,2)
plot(xvec/pi,(Vsl1+Vsl2)*0.5,'b-',xvec/pi,0.5*Vsl1,'r--',xvec/pi,0.5*Vsl2,'r--',xvec/pi,rhox,'k-')
xlabel('x/\pi'); ylabel('V(x) & \rho(x)');
hold on
title(['N=',num2str(N),'   \eta = ',num2str(Eta)]);











