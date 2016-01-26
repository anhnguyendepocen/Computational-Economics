function [x,fx,ef,iter] = ps244(f,x,I)
[value] = f(x); fx=value(:,1);
n = length(x); Jac = eye(n,n);
for j = 1:I(1,3) %broyden method
    xp = x - Jac\fx;
    [valuep]= f(xp); fxp=valuep(:,1);
    xd = xp-x; fxd = fxp - fx;
    Jacp = Jac + (fxd-Jac*xd)*xd' / (xd'*xd);
    D = (norm(x-xp) <= I(1,1)*(1+norm(xp)) && norm(fx) <=I(1,2));
    if D == 1; 
        break;
    else
        x = xp; fx = fxp; Jac = Jacp; 
    end
end
ef = 0; if D == 1; ef = 1; end
iter = j;