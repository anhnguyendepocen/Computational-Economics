function [x,fx,ef,iter] = ps243(f,x,I)
for j = 1:I(1,3)%Newton Method
    value=f(x);
    fx=value(:,1);dfx=value(:,2:3);
    xp = x - dfx\fx;
    D = (norm(x-xp) <= I(1,1)*(1+norm(xp)) && norm(fx) <=I(1,2));
    if D == 1; 
        break;
    else
        x = xp;
    end
end
ef = 0; if D == 1; ef = 1; end
iter = j;