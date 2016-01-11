function y = cournot(x)

global n lambda zeta

y = n - 1/lambda;

for i=1:n
    y = y - 1/(x^(-1-1/lambda)/lambda/zeta(i) + 1);
end;

end
