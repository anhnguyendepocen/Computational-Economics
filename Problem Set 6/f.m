function euler_resid = f(k_tild)

global p alpha beta gamma sigma lam gamma_tild k g;

phi1=(((1-alpha)*((1-alpha)/(1+lam))^(alpha/(1-alpha)))^(1-sigma))/(1-sigma);
phi_tild=(((1+lam)/(1-alpha))^(alpha*(1-sigma)))/((1-alpha)^(1-sigma));
A=(1-alpha)/(1+lam);
s= k_tild/(A*k(p)^alpha);
g_hat=interp1(k,g,k_tild,'spline');
s_prime=g_hat/(A*k_tild.^alpha);
m=alpha*(1-sigma)/(1-alpha);
Vs=phi1*(-(1-sigma)*(1-s)^(-sigma))*s^m+(((1-s)^(1-sigma)+beta*phi_tild*gamma_tild)*m*s^(m-1));
Vs_prime=phi1*(-(1-sigma)*(1-s_prime)^(-sigma))*s_prime^m+(((1-s_prime)^(1-sigma)+beta*phi_tild*gamma_tild)*m*s_prime^(m-1));
euler_resid =Vs_prime - (1/gamma+1)*Vs;
end