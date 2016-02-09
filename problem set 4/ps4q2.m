clc; close all; clear all;
%% Problem Set 4, Question 2
%Set parameters
w=(0.5:0.5:50)';rmin=-0.08;rhigh=0.12;gemma=2;p=0.5;%for rmin
%linear solution
c_l=0.5*w*(1+p*rmin+(1-p)*rhigh);
%c_che=ones(1,length(w));c00=1;
%CRRA
%there is some problem to solve in this way, however fine with the second
% syms c y
% y=c.^(-gemma)-p.*(((1+rmin).*w-c).^(-gemma))-(1-p).*(((1+rhigh).*w-c).^(-gemma));
% f=matlabFunction(y);
% c_che=fsolve(f,c0);
%c_che=fsolve(@(c)c.^(-gemma)-p.*(((1+rmin).*w-c).^(-gemma))-(1-p).*(((1+rhigh).*w-c).^(-gemma)),c0);
%it takes some while to solve that,around 30sec but the answer is terrible
%plot(w,c_l,w,c_che)

% To solve consumption at very wealth point instead of vectorize
% But the answer doesn't improve with both guesses
% for i=1:1:length(w)
% c_che(i)=fsolve(@(c)c.^(-gemma)-p.*(((1+rmin).*w(i)-c).^(-gemma))-(1-p).*(((1+rhigh).*w(i)-c).^(-gemma)),c00);
% end
% plot(w,c_l,w,c_che)
% 
% figure
% for i=1:1:length(w)
% c_che(i)=fsolve(@(c)c.^(-gemma)-p.*(((1+rmin).*w(i)-c).^(-gemma))-(1-p).*(((1+rhigh).*w(i)-c).^(-gemma)),c_l(i));
% end
% plot(w,c_l,w,c_che)
% And then comes the meaning of interpolation, as we don't know the closed
% form answer for the question we are going to solve

%Chebychev intepolation based on M&F collation method
n=44;a=w(1);b=w(length(w));%critical change when n=44
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

%maximum percentage error of the deviation.
mped=max(abs(splot./c2));%actually it is a one vector

% %Chebychev intepolation based on M&F collation method
% m=44;a=w(1);b=w(length(w));%critical change when n=44
% fspace_m = fundefn('cheb',m,a,b);
% w_ch_m = funnode(fspace_m);
% coeff_m = 0.1*(1:1:length(w_ch_m))';
% %coeff = ones(length(w_ch),1);
% coeff_m = broyden('resid',coeff_m,w_ch_m,fspace_m,gemma,p,rmin,rhigh);
% splot_m = funeval(coeff_m,fspace_m,w);
% figure
% subplot(1,2,1)
% plot(w,splot_m);%plot the risidual function
% xlim([0.05,45])
% c_m=funeval(coeff_m,fspace_m,w_ch_m);
% subplot(1,2,2)
% plot(w_ch_m,c_m)
% xlim([0.05,45])
% xlabel('Wealth')
% ylabel('Consumption')

