clc
close all; clear all;

global n lambda zeta q aq
syms x y
q = zeros(1,n); aq = 0;

n = input('The number of firms:');
lambda = 1.6;
zeta = linspace(0.6, 0.8, n);

y = n - 1/lambda;
for i=1:n
    y = y - 1/(x^(-1-1/lambda)/lambda/zeta(i) + 1);
end;
cournot = matlabFunction(y);

aq = fzero(cournot, sqrt(n)*2); % Aggregate production

if (aq<=0) disp('Non-positive aggregate production!'); return; end;

for i=1:n
    q(i) = aq^(-1/lambda)/(zeta(i)+aq^(-1-1/lambda)/lambda); % Each firm's production
end

disp(['The equilibrium allocations are: ', num2str(q,' %.4f ')]);
