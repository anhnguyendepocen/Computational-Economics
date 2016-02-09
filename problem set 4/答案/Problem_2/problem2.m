% problem 2
function problem2

close all

global gam

rf = 0.02;
er = 0.1;
r = [rf-er;rf+er];
p = 0.1;
gam=2;     % coefficient of risk aversion
c=1-gam;

% for exact solution:
nlarge=501;

er = r-rf;
prob = [p; 1-p];
RF = 1+rf;

w0=0.5;
w1=50;
n=15;

fspace=fundefn('cheb',n,w0,w1);
wgrid=funnode(fspace);

% linear function as starting values:
clin=prob'*(1+r)/2*wgrid;
B=funbas(fspace,wgrid);
coefflin=B\clin;

% coefficients:
coeff=broydn(@eulerres,coefflin,1e-6,0,0,fspace,wgrid,prob,r);
disp('old coefficients: ');
coefflin
disp('new coefficients: ');
coeff
cgrid=B*coeff;

% accuaracy: Euler equation errors:
wlarge=nodeunif(nlarge,w0,w1);
res=feval(@eulerres,coeff,fspace,wlarge,prob,r);
disp(['maximum euler equation errors: ', num2str(max(abs(res)))]);

figure;
plot(wgrid, clin, 'r--', wgrid, cgrid, 'b-');
legend('linear approx');
disp(['maximum percentage deviation: ', num2str(max(abs(cgrid-clin)./clin))]); 


% -------------------------------------------------------------------------
function res=eulerres(coeff,fspace,w,prob,r)

global gam

cons=funeval(coeff,fspace,w);
Eucp1=0;
for i=1:length(prob),
    Eucp1=Eucp1+prob(i)*(w.*(1+r(i))-cons).^-gam;
end;
res=1-Eucp1.^(-1/gam)./cons;
% -------------------------------------------------------------------------

