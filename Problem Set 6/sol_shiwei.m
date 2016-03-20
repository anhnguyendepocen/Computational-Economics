clc; clear all; close all;
%Howard Improvement
%% Initialization
% parameters setting
%alpha, beta, lamda, sigma, gamma,tilda
a=0.33; b=0.7; l=0.5; s=2; g=0.7; t=1;
phi_1=(((1-a)*((1-a)/(1+l))^(a/(1-a)))^(1-s))/(1-s);
phi_2=(((1-a)/(1+l))^((a^2*(1-s))/(1-a)))/(1-s);
ind=a*(1-s)/(1-a);
%value function for saving rate 
%v(sr)=(phi_1*(1-sr)^(1-s)+b*phi_2*t)*s^ind

%value function for k
%w(k)=max(s){v(k,sr)+g*w(k')}
%k'=(1-a)/(1+l)*sr*k^a
%v(k,sr)=(((1-sr)*(1-a)*k^a)^(1-s)+b*(k^a*sr)^(a*(1-s))*t)/(1-s)

%Other setting
k0=(1:1:100)/20;tol=1.0e-12;imax=100;	
n=length(k0);k=k0;
w0=ones(1,n);w=zeros(1,n);index=zeros(1,n);wnew=w0;
%kgrid1=repmat(fliplr(k0),n,1)';kgrid2=repmat(k0,n,1)
%% Basic value iteration
j=1;
while max(abs(wnew-w))>tol&&j<imax
w=wnew;
for m=1:n
sr = k*(1+l)/(1-a)/k(m)^a;%1*n dimension
v=value(k(m),sr,a,b,s,t,n);%1*n
wtemp=v + g*w;%1*n
% In the next line, index is position of max element
[wnew(m),index(m)]=max(wtemp);
end
j=j+1;
end
kstar=k(index);
srstar = kstar*(1+l)/(1-a)./(k(m)^a);%1*n dimension
plot(k,srstar,k,kstar)
xlabel ('k'); ylabel ('ratio ');
legend('saving rate','k_(prime)','Location','northwest')
for m=1:n;
    vstar=value(kstar(m),srstar,a,b,s,t,n);
end
figure
plot(k,w,k,vstar)%seems a problem in vstar, which is value function in the paper
xlabel ('k'); ylabel ('utility');
legend('social welfare','value','Location','northwest')

%% Howard improvement w0=ones not good one
% kh=0.2*w0;knew=k0;indexh=index;j=1;
% syms kstarh 
% while max(abs(knew-kh))>tol||j>imax
%     kh=knew;
% for m=1:n
% sr = kstarh*(1+l)/(1-a)/k(m)^a;%1*n dimension
% v=valueh(k(m),sr,a,b,s,t,n);%1*n
% wtemp1=1/(1-g)*v;%1*n
% wtemp2=v+g*wtemp1;
% y= matlabFunction(-wtemp2);
% %knew(m)=fminsearch(y,kh(m));
% knew(m)=fmincon(y,kh(m),-1,0);
% end
% j=j+1;
% end
% sr= knew*(1+l)/(1-a)./(k(m)^a);
% figure
% plot(k,knew,k,sr)
% xlabel ('k'); ylabel ('ratio ');
% legend('k_(prime)','saving rate','Location','northwest')

%% policy Iteration
% w0=ones 
j=1;knew=w0;k=0*w0;indexh=index;%knew is intial guess
while max(abs(knew-k))>tol&&j<imax
    k=knew;%updating policy guess
    sr=k*(1+l)/(1-a)./k0.^a; %updating policy guess
    v=value(k0,sr,a,b,s,t,n);%updating value
    wtemp1=1/(1-g)*v;%updating value
for m=1:n; %find knew 
    sr=k0*(1+l)/(1-a)./k0(m)^a;
    v=value(k0,sr,a,b,s,t,n);
     wtemp2=v+g*wtemp1;
[wnew(m),indexh(m)]=max(wtemp2);
knew(m)=k0(indexh(m));
end
j=j+1;
end
kstarh=k0(indexh);
srstarh = kstarh*(1+l)/(1-a)./(k0(m)^a);%1*n dimension
figure
plot(k0,srstarh,k0,kstarh)
xlabel ('k'); ylabel ('ratio ');
legend('saving rate','k_(prime)','Location','northwest')
vstarh=value(k0,srstarh,a,b,s,t,n);
figure
plot(k0,wnew,k0,vstarh)%seems a problem in vstar, which is value function in the paper
xlabel ('k'); ylabel ('utility');
legend('social welfare','value','Location','northwest')