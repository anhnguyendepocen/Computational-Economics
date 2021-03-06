clc
close all; clear all;

syms x y;

y = x^3 + 4 - 1/x;
g1 = matlabFunction(y);
[x1, f1] = bisection(g1, [0.1; 1], 1e-5, 1e5);
disp(['The root of problem 2.1 is: ', num2str(x1,' %.4f ')]);

y = -exp(-x) + exp(-x^2);
g2 = matlabFunction(y);
[x2, f2] = bisection(g2, [-0.5; 0.5], 1e-8, 1e8);
disp(['The root of problem 2.2 is: ', num2str(x2,' %.4f ')]);

a = 3; b = 0.5; c = 1; d = 1; phi=0.5;
y = b*x + d*x^phi - (a-c);
ds = matlabFunction(y);
[q1,qf1] = bisection(ds, [0; 2], 1e-8, 1e8);
p1 = a - b*q1;
disp(['The solution of problem 3 with the bisection method is: ', ' q = ', num2str(q1,' %.4f, '), ' p = ', num2str(p1,' %.4f ')]);
q2 = fzero(ds, 1);
p2 = a - b*q2;
disp(['The solution of problem 3 with fzero command is: ', ' q = ', num2str(q2,' %.4f, '), ' p = ', num2str(p2,' %.4f ')]);
