% problem3

clear
close all
clc

rf = 0.02;
er = 0.06;
r = [rf-er;rf+er];
p = 0.1;

n=15;

prob = [p; 1-p];

gam0=1;
gam1=50;
addpath(genpath('CEdemos'));
addpath(genpath('CEtools'));

fspace=fundefn('cheb',n,gam0,gam1);
gamgrid=funnode(fspace);

alph=zeros(n,1);
alph0=alph(1);
for i=n:-1:1,
    alph(i)=broyden(@alphres,alph0,gamgrid(i),prob,r,rf);
    alph0=alph(i);
    if alph0>1,
        gam=gamgrid(i);
        alph(1:i)=1;
        
        break;
    end;
end;
B=funbas(fspace,gamgrid);
coeff=B\alph;

figure;
plot(gamgrid,alph,'r-');

% smarter way:
gam0new=fzero(@resgam0,[gam0; gam1],[],prob,r,rf);

fspace=fundefn('cheb',n,gam0new,gam1);
gamgrid=funnode(fspace);

% coefficients:
alphgrid=broyden(@alphres,ones(n,1),gamgrid,prob,r,rf);

gamgrid=[gam0;gamgrid];
alphgrid=[1;alphgrid];

hold on;
plot(gamgrid,alphgrid,'y--');

