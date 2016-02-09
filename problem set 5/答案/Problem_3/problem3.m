% problem3
function problem3

clc
close all

n=11;       % number of nodes
sig=0.25;   % s.d.
gam=1.0;    % risk aversion
Rf=1.02;
R=1.1;
epsi=1.0e-6;

% approximate log-normally distributed random variable
sig2=sig^2;
mu=-sig2/2;
[x,w]=qnwnorm(n,mu,sig2);
x=exp(x);

% test
ex=w'*x;
dist=ex-1.0;
if (abs(dist)>epsi)
    disp('increase number of nodes');
    return;
end;

% form stochastic (gross) return
R=R*x;

% specify functions:
utilfun=@(gam,c) f_util(gam,c);
distfun=@(gam) f_dist(utilfun,gam,Rf,R,w);

% TEST
% % utility for risk-free investment
% util_rf=feval(utilfun,gam,Rf); 
% 
% % utility for risky investment
% util_r=feval(utilfun,gam,R);
% util_r=w'*util_r;

% test: plot distance function
gamvec=[0.1:0.1:4]';
ng=length(gamvec);
dist=zeros(ng,1);
for gc=1:ng,
    dist(gc)=feval(distfun,gamvec(gc));
end;
plot(gamvec,dist);

% evaluate distance at baseline parameter
dist=feval(distfun,gam);
if (abs(dist)>epsi),
    gam=fzero(distfun,gam);
    disp(['solution for gamma is: ', num2str(gam)]);
else
    disp('no numerical solution is required');
end;


end

function util=f_util(gam,c)

epsi=1.0e-04;

if (abs(gam-1.0)<epsi),
    util=log(c);
else
    temp=1.0-gam;
    util=1.0./temp*(c.^temp-1.0);
end

end

function dist=f_dist(utilfun,gam,Rf,R,w)

f1=feval(utilfun,gam,Rf);
f2=feval(utilfun,gam,R);
f2=w'*f2;

dist=f1-f2;

end