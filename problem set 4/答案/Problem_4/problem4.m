%% Masterfile for Problem 4
close all; clear all 
clc

%pathRoot = fileparts(which(mfilename));
%addpath(genpath(pathRoot)); 
cputime = 0; tic;


%% Step 01: Parameterization
globalvars;
alp = 0.36;     % capital share in output
bet = 0.96; 	% time preference rate
delt = 0.06;    % depreciation rate



%% Step 02: Steady State Capital Stock and Grid
KSS = ( (1/bet+delt-1)/alp )^(1/(alp-1));
Kmin = 0.5*KSS; Kmax = 1.5*KSS;

m = 9; n = m-1;
normgrid = -cos( (2*[1:m]-1)/(2*m) * pi )'; 
[Tp] = poly1D_cheb(n,normgrid);



%% Step 03: Euler Equation and Root Finding;
% generate initial guess for coefficents: suppose K = KP, then, C = (r-delt)*K+w
K = (normgrid*(Kmax-Kmin)+Kmin+Kmax)/2;     % real capital grid in t
r = alp*K.^(alp-1);                         % interest rate in t
w = (1-alp)*K.^(alp);                       % wage rate in t
C0 = (r-delt).*K+w;                         % consumption when K = Kp    
coeff0 = Tp\C0;                             % coefficents for C0
if n > 20;
    load('coeffini');
    coeff0(1:21,1) = coeff11;
end
    
% root finding algorithm
options_fsolve = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12, 'MaxIter', 3000, 'MaxFunEvals', 3000);
f = @(x) eulereq(x,normgrid,Tp);
[x01,fv01,ef01,op01,jac01] = fsolve(f,coeff0,options_fsolve);
if ef01 ~= 1
    error('rootfinding algorithm failed, exitflag = %.1d\n',ef01);4
end
coeff = x01;
if n == 20;
    coeff11 = coeff; 
    save('coeffini.mat','coeff11');
end


%% Step 04: Plot
plotgrid=linspace(-1, 1, 100)';
[Tp2] = poly1D_cheb(n,plotgrid);
plotgrid=(plotgrid +1)*(Kmax-Kmin)/2+Kmin;
Cp = Tp2*coeff;
figure;
plot(plotgrid,Cp,'x');

%% Step 05: Simulate capital
nt=100;
Ct=zeros(nt,1);
Kt=zeros(nt,1);
Kt(1)=Kmin;

for t=1:nt-1
    rt = alp*Kt(t).^(alp-1);                         % interest rate in t
    wt = (1-alp)*Kt(t).^(alp);                       % wage rate in t
    K_scaled=2*(Kt(t)-Kmin)/(Kmax-Kmin)-1;
    Ct(t)= poly1D_cheb(n,K_scaled)*coeff;
    Kt(t+1)=(1+rt-delt)*Kt(t)- Ct(t) +wt;
end
Ct(nt)=poly1D_cheb(n,Kt(nt))*coeff;
figure;
plot([1:100],Kt,'x');

%% Step 06: Compute Euler Error on Random Sample of Points
randgrid = 2*(sortrows(rand(1000,1))-0.5);
[Tp] = poly1D_cheb(n,randgrid);
[err,C,Cp] = eulereq(coeff,randgrid,Tp);

% relative Euler error
releerr = (Cp-C)./C; 
fprintf('mean relative Euler error in percent:  %8.4f \n', 100*mean(abs(releerr)))
fprintf('max relative Euler error in percent:   %8.4f \n', 100*max(abs(releerr)))


