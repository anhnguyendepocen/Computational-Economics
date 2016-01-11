function[x, fx] = bisection(f, x, crit, niter)

xl = x(1,1); fl = f(xl);
xh = x(2,1); fh = f(xh);

if (fh*fl > 0)
    disp('initial [xl, xh] do not bracket a root');
    return
end

iter = 0;

while (xh-xl)>crit*(1+abs(xl)+abs(xh)) && iter<=niter
    iter = iter + 1;
    xm = (xl+xh)/2; fm = f(xm);
    if (fm==0) break; end;
    if fl*fm < 0
        xh = xm; fh = fm;
    else
        xl = xm; fl = fm;
    end
end

x = xm; fx = fm;
