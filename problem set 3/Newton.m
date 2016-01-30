function [x,fx,ef] = Newton(f,x,cc,t)
tole = cc(1,1); told = cc(2,1); maxiter = cc(3,1);
ef = 0;x0=x;
for j = 1:maxiter
[fx,dfx,ddfx] = f(x);
xp = x - ddfx\dfx';
if norm(x-xp) <= tole*(1+norm(xp))
ef = 2; break; % converged in first criterion
else
x = xp;
end
end
if norm(dfx) <= tole*(1+norm(fx))
if ef == 2; ef=1; % converged in both criteria
else ef=3; % spurious convergence
end
end
sol=[ef,x,fx];
fprintf('for starting value %4.2f of this funtion %s',x0,t)
disp('')
fprintf(' the solution indicator is %1.0f,solution x is %8.3f, function value fx is %8.3f', sol)
disp('')

