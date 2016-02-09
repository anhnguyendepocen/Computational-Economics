clc; close all; clear all;
%Part two with different parameters
%1.increas gemma, was 2;2.decrease p, was 0.5;3.enlarge the gap of rmin&rhigh;was
%-0.08&0.12
w=(0.5:0.5:50)';rmin=-0.08;rhigh=0.12;gemma=10;p=0.5;%for rmin
%linear solution
c_l=0.5*w*(1+p*rmin+(1-p)*rhigh);
%Chebychev intepolation based on M&F collation method
n=15;a=w(1);b=w(length(w));%critical change when n=44
fspace = fundefn('cheb',n,a,b);
w_ch = funnode(fspace);
coeff = 0.1*(1:1:length(w_ch))';
%coeff = ones(length(w_ch),1);
coeff = broyden('resid',coeff,w_ch,fspace,gemma,p,rmin,rhigh);
splot = funeval(coeff,fspace,w);
plot(w,splot);%plot the risidual function
c1 = funeval(coeff,fspace,w_ch);
c2 = funeval(coeff,fspace,w);
figure
plot(w,c_l,'--',w,c2,w_ch,c1,'-o')
legend('linear function','chebychev on whole x','chebychev on nodes','Location','NorthWest')
xlabel('Wealth')
ylabel('Consumption')