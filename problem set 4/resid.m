function r = resid(coeff,w_ch,fspace,gemma,p,rmin,rhigh)
c = funeval(coeff,fspace,w_ch);
r = -c+(p.*(w_ch.*(1+rmin)-c).^(-gemma)+(1-p).*(w_ch.*(1+rhigh)-c).^(-gemma)).^(-1/gemma);
%very sensitive to the form of the residual, for example c or c^-gemma on
%the LHS