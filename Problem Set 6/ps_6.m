clear
close all
clc


%% Baseline Calibration
alpha=0.33;                 % output elasticity of capital
beta=0.7;                  % discount factor of households
gamma=0.7;                  % discount factor of government
sigma=2.0;                  % coefficient of relative risk aversion
lam=0.5;                    % second period productivity
sig_eta=0.65;% standard deviation of idiosyncratic income shocks
nk =50;
k=linspace(0.3,0.45,nk)';    % number of grid points in capital grid
% ns=10;
% s=linspace(0.01,0.45,ns)';    % number of grid points in saving rate grid (pick maximum)
tolf=1.0e-12;               % tolerance for function distance
g=ones(nk,1);
gnew=linspace(0.3,0.45,nk)';
gamma_tild=1;

%% Iteration on the first-order conditions
options_fsolve = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', 3000, 'MaxFunEvals', 3000);
while max(abs(gnew-g))>tolf
 g=gnew;
for p=1:nk
 gnew(p) = fsolve(@f,g(p));
end
end
poli_k=g;
plot(k,poli_k)


































