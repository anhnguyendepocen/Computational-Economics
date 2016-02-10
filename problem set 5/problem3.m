% problem3


clc
close all
addpath(genpath('CEdemos'));
addpath(genpath('CEtools'));
n=11;       % number of nodes
sig=0.25;   % s.d.
gam=1.5;    % risk aversion
y1=1.02;
y2=1.06;
y2=1.1
epsi=1.0e-6;

% approximate log-normally distributed random variable
sig2=sig^2;
mu=-sig2/2;
[x,w]=qnwnorm(n,mu,sig2);
elnx=w'*x;
x=exp(x);

% test
ex=w'*x;
dist=ex-1.0;
if (abs(dist)>epsi)
    disp('increase number of nodes');
    return;
end;

% form stochastic (gross) return
R=y2*x;
% expected return of two projects
ey1=y1
ey2=ex*y2
% specify functions:
utilfun=@(gam,c) f_util(gam,c);
distfun=@(gam) f_dist(utilfun,gam,y1,R,w);

% TEST
% utility for risk-free investment
util_rf=feval(utilfun,gam,y1); 

% utility for risky investment
util_r=feval(utilfun,gam,y2);
util_r=w'*util_r;

% test: plot distance function
gamvec=[0.1:0.1:4]';
ng=length(gamvec);
dist=zeros(ng,1);
for gc=1:ng,
    dist(gc)=feval(distfun,gamvec(gc));
end;
plot(gamvec,dist);
xlabel('gamma')
ylabel('Difference between two utilities')
title('Distance Function when y2=1.06')

% evaluate distance at baseline parameter
dist=feval(distfun,gam);
if (abs(dist)>epsi),
    gam=fzero(distfun,gam);
    disp(['solution for gamma is: ', num2str(gam)]);
else
    disp('no numerical solution is required');
end;



