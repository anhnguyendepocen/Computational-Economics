function [err,C,Cp] = eulereq(coeff,normgrid,Tp)


globalvars;


n = length(coeff)-1;
K = (normgrid*(Kmax-Kmin)+Kmin+Kmax)/2;     % real capital grid in t

% interets and wage rate
r = alp*K.^(alp-1);                         % interest rate in t
w = (1-alp)*K.^(alp);                       % wage rate in t

% consumption and saving in period t
C = Tp*coeff;                               % consumption in t
Kp = (1+r-delt).*K + w - C;                 % capital stock in t+1

% consumption and interest rate in period t+1
normgridp = (2*Kp-Kmin-Kmax)/(Kmax-Kmin);   % normalized capital grid in t+1
[Tpp] = poly1D_cheb(n,normgridp);
Cp = Tpp*coeff;
rp = alp*Kp.^(alp-1);

% Euler equation
err = Cp-bet*(1+rp-delt).*C;
