% problem3
function problem3

close all


rf = 0.02;
er = 0.1;
r = [rf-er;rf+er];
p = 0.1;

n=11;
nlarge=501;

er = r-rf;
prob = [p; 1-p];

gam0=1;
gam1=50;

fspace=fundefn('cheb',n,gam0,gam1);
gamgrid=funnode(fspace);

alph=zeros(n,1);
alph0=alph(1);
for i=n:-1:1,
    alph(i)=broydn(@alphres,alph0,1e-6,0,0,gamgrid(i),prob,r,rf);
    alph0=alph(i);
    if alph0>1,
        alph(1:i)=1;
        break;
    end;
end;
B=funbas(fspace,gamgrid);
coeff=B\alph;

figure;
plot(gamgrid,alph,'b-');

% smarter way:
gam0new=fzero(@resgam0,[gam0; gam1],[],prob,r,rf);

fspace=fundefn('spli',n,gam0new,gam1);
gamgrid=funnode(fspace);

% coefficients:
alphgrid=broydn(@alphres,ones(n,1),1e-6,0,0,gamgrid,prob,r,rf);

gamgrid=[gam0;gamgrid];
alphgrid=[1;alphgrid];

hold on;
plot(gamgrid,alphgrid,'r--');

% -------------------------------------------------------------------------
function res=alphres(alph,gam,prob,r,rf)

res=0;
for i=1:length(prob),
    res=res+prob(i)*(1+rf+alph.*(r(i)-rf)).^(-gam).*(r(i)-rf);
end;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function res=resgam0(gam,prob,r,rf);

res=prob'*((1+r).^(-gam).*(r-rf));