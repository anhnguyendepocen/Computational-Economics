% Problem 4.2 and 4.3

syms a;
rf = 0.02; rlr = -0.08-rf; rhr = 0.12-rf;
p = 0.1; phi = -3;

f = -(p*(1 + rf + a*rlr)^phi + (1 - p)*(1 + rf + a*rhr)^phi)/phi;
df = diff(f,a);
mf = matlabFunction(f);
mdf = matlabFunction(df);
astar1 = fzero(mdf,1);
display(['The unconstrained solution is ', num2str(astar1)]);

figure;
set(gcf,'color','w');
ezplot(-f,[astar1-1,astar1+1]);
xlabel('value of alpha');
ylabel('value of the objective (w0=1)');
title('');
print('plot42','-dpng')

astar2 = fminbnd(mf,0,1);
astar3 = fmincon(mf,0.5,[],[],[],[],0,1);
display(['The function fminbnd and fmincon give the constrained solution ', num2str(astar2), ' and ', num2str(astar3)]);