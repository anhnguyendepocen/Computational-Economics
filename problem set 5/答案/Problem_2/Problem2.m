% problem2
function Problem2

close all

global gam

disp('solution to this problem calls function kronmex');
disp('kronmex calls a function kronc, coded in C');
disp('to use kronc, you will have to first compile it');
disp('to do this, type mex kronc.c');
inpt=input('have you already mexed kronc? If yes, press 1: ');
if (inpt~=1),
    disp('execution terminated');
    return
end;


% settings
gam=1;
w0 = 100;
wmin = 20;
rf = 0.02;
r = [0.04; 0.06];
rho = 0.5;
Std = [0.1; 0.2];
Cov = rho*Std(1)*Std(2);
VarCov = diag(Std.^2);
VarCov(1,2)=Cov;
VarCov(2,1)=Cov;
if any(eig(VarCov)<0),
    error('variance-covariance matrix must be p.d.');
end;
n = 3;

[rnodes, rwghts] = qnwnorm([n, n], r, VarCov);

alpha0=ones(length(r),1);
res = alphres(alpha0,w0,wmin,rf,rnodes,rwghts);

% unconstrained problem:
disp(' ');
disp('unconstrained problem:');
for wmin=0:10:50,
    alph=broyden(@alphres,alpha0,w0,wmin,rf,rnodes,rwghts);
    alpha0=alph;
    disp(['pf-shares for minimum wealth level ', num2str(wmin), ' are: ', num2str(alph')]); 
end;

% constrained problem:
optset('ncpsolve','type','smooth');
n=length(r);
alpha0=ones(length(r),1);
disp(' ');
disp('constrained problem:');
for wmin=0:10:50,
    alph=ncpsolve(@ncpalphres,zeros(n,1),ones(n,1),alpha0,w0,wmin,rf,rnodes,rwghts);
    alpha0=alph;
    disp(['pf-shares for minimum wealth level ', num2str(wmin), ' are: ', num2str(alph')]); 
end;


% -------------------------------------------------------------------------
function [fval,fjac]=ncpalphres(alph,w0,wmin,rf,rnodes,rwghts);

fval = feval(@alphres,alph,w0,wmin,rf,rnodes,rwghts);
fjac = fdjac(@alphres,alph,w0,wmin,rf,rnodes,rwghts);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function res=alphres(alph,w0,wmin,rf,rnodes,rwghts);

global gam;

n=length(alph);
res=zeros(n,1);
for i=1:n,
    res(i) = rwghts'*((w0.*(1+rf+(rnodes-rf)*alph)-wmin).^-gam.*w0.*(rnodes(:,i)-rf));
end;
% -------------------------------------------------------------------------


