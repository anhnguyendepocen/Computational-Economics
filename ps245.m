function [x,fx,ef,iter] = ps245(f,x,I)
% inverse broyden
[value] = f(x); fx=value(:,1);
n = length(x); B = eye(n,n);
for j = 1:I(1,3)
    xp = x - B*fx;
    [valuep]= f(xp); fxp=valuep(:,1);
    xd = xp-x; fxd = fxp - fx;
    Bp = B + (xd-B*fxd)*xd'*B / (xd'*B*fxd);
    D = (norm(x-xp) <= I(1,1)*(1+norm(xp)) && norm(fx) <=I(1,2));
    if D == 1; 
        break;
    else
        x = xp; fx = fxp; B = Bp; 
    end
end
ef = 0; if D == 1; ef = 1; end
iter = j;