clc; clear all; close all;

%% Initialization
syms x;
y = [x^4, x^6, 1/(1+x^2)]';
my = matlabFunction(y);

%% Guassian-Hermite Quadrature Integration
n= [2, 3, 4, 5, 7];
for i = 1:length(n)
    [x, w] = qnwnorm(n(i), 0, 1); % by Miranda and Fackler (2002)
    int = sum(repmat(w',3,1).*my(x'),2);
    disp(sprintf('When n = %d, GH integrations are %4.3f, %4.3f, %4.3f', n(i), int));
end

%% Monte Carlo Integration
n = [1e2, 1e3, 1e4, 5*1e4];
for i =1:length(n)
    r = randn(1,n(i));
    int = sum(my(r),2)./n(i);
    disp(sprintf('When n = %d, MC integrations are %4.3f, %4.3f, %4.3f', n(i), int));
end