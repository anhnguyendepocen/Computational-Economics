close all; 
clear all 
clc;
%% PS3 q1 Newton Method
cc=[eps,eps,100]';x0=2;

syms x y ;
y = 2*x^3-x^2-3*x+2;
dy=diff(y);
ddy=diff(dy);
f=matlabFunction(y,dy,ddy);
txt='y = 2*x^3-x^2-3*x+2';
fsol=Newton(f,x0,cc,txt);

disp(' ')

y=-x*exp(-x);
dy=diff(y);
ddy=diff(dy);
g=matlabFunction(y,dy,ddy);
txt='y = -x*exp(-x)';
gsol=Newton(g,x0,cc,txt);

% solution with Matlab solver:
y=2*x^3-x^2-3*x+2;
f1= matlabFunction(y);
xmnew_f = fminunc(f1,x0);

y=-x*exp(-x);
g1= matlabFunction(y);
xmnew_g = fminunc(g1,x0);
