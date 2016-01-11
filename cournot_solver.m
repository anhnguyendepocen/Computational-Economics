clc
close all; clear all;

global n lambda zeta q aq sum

n = input('The number of firms:');
lambda = 1.6;
zeta = linspace(0.6, 0.8, n);

q = zeros(1,n); aq = 0; sum = 0;

aq = fzero(@cournot, sqrt(n)*2);

if (aq<=0) disp('Non-positive aggregate production!'); return; end;

for i=1:n
    q(i) = aq^(-1/lambda)/(zeta(i)+aq^(-1-1/lambda)/lambda);
    sum = sum + q(i);
end

disp(['The equilibrium allocations are: ', num2str(q,' %.4f ')]);
