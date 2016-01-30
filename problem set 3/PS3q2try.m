clc; close all; clear all;
%%Multi equation 
%set parameters and variables
beta=0.95;r=0.02;theta=1;T=10;max=1000;tol=1e-30;
c=ones(1,T);a0=c;w=zeros(1,T);w(1)=10;t=(1:1:T);a0(1)=0;
%a(1)=0,a(t+1)=0,a can't be set as all zeros in the first time

%(a(T-1)*(1+r)+w(T-1)-a(T))^(-theta)-beta*(1+r)*(a(T)*(1+r)+w(T))^(-theta)=0;
%FOC for T-1
%further with respect to a(T)
%theta*(a(T-1)*(1+r)+w(T-1)-a(T))^(-theta-1)+theta*beta*(1+r)^2*(a(T)*(1+r)+w(T))^(-theta-1)=0
%for i=T-2:-1:1 (a(i)*(1+r)+w(i)-a(i+1))^(-theta)-beta*(1+r)*(a(i+1)*(1+r)+w(i+1)-a(i+2))^(-theta)=0
%rest of the period
%further with respect to a(i+1)
%theta*(a(i)*(1+r)+w(i)-a(i+1))^(-theta-1)+theta*beta*(1+r)^2*(a(i+1)*(1+r)+w(i+1)-a(i+2))^(-theta-1)=0

%Linear Gauss-Jacobi algorithm:
%xi_k+1=xi_k-gi(xi_k)/gi'(xi_k)£¬no updating further
a=a0;a1=a0;
tic;
for j=1:max
    for i=T-1:-1:2
        a1(T)=a(T)-((a(T-1)*(1+r)+w(T-1)-a(T))^(-theta)-beta*(1+r)*(a(T)*(1+r)+w(T))^(-theta))/(theta*(a(T-1)*(1+r)+w(T-1)-a(T))^(-theta-1)+theta*beta*(1+r)^2*(a(T)*(1+r)+w(T))^(-theta-1));
        a1(i)=a(i)-((a(i-1)*(1+r)+w(i-1)-a(i))^(-theta)-beta*(1+r)*(a(i)*(1+r)+w(i)-a(i+1))^(-theta))/(theta*(a(i-1)*(1+r)+w(i-1)-a(i))^(-theta-1)+theta*beta*(1+r)^2*(a(i)*(1+r)+w(i)-a(i+1))^(-theta-1));
        %a1(1)=a(1)-((-a(1))^(-theta)-beta*(1+r)*(w(1)-a(1))^(-theta))/(theta*(-a(1))^(-theta-1)+theta*beta*(1+r)^2*(a(1)*(1+r)+w(1)-a(2))^(-theta-1));
    end
    dist=abs(a1-a);
    if (dist<tol), 
        break
    else
        a=a1;
    end
end
t1=toc;a1
fprintf('the time for Gauss-Jacobi algorithm is %8.6f second',t1)
disp(' ')
%consumption setting;
c1=c;
for i=1:1:T-1;
 c1(i)=a1(i)*(1+r)+w(i)-a1(i+1);
 end
 c1(T)=a1(T)*(1+r)+w(T);
 c1
plot(t,c1,'r',t,w,'-*',t,a1,'--o')
title('PS3,Q2,Consumption Saving Problem Linear Gauss-Jacobi algorithm')
xlabel('time period')
xlim([1 T])
ylabel('Amount')
legend('c','y','a','Location','NorthEast')

%Linear Gauss-Siedel algorithm:
%xi_k+1=xi_k-gi(xi_k)/gi'(xi_k)£¬no updating further
a=a0;a2=a0;
tic
for j=1:max
 
    for i=T-1:-1:2
        a(T)=a(T)-((a(T-1)*(1+r)+w(T-1)-a(T))^(-theta)-beta*(1+r)*(a(T)*(1+r)+w(T))^(-theta))/(theta*(a(T-1)*(1+r)+w(T-1)-a(T))^(-theta-1)+theta*beta*(1+r)^2*(a(T)*(1+r)+w(T))^(-theta-1));
        a(i)=a(i)-((a(i-1)*(1+r)+w(i-1)-a(i))^(-theta)-beta*(1+r)*(a(i)*(1+r)+w(i)-a(i+1))^(-theta))/(theta*(a(i-1)*(1+r)+w(i-1)-a(i))^(-theta-1)+theta*beta*(1+r)^2*(a(i)*(1+r)+w(i)-a(i+1))^(-theta-1));
        %a1(1)=a(1)-((-a(1))^(-theta)-beta*(1+r)*(w(1)-a(1))^(-theta))/(theta*(-a(1))^(-theta-1)+theta*beta*(1+r)^2*(a(1)*(1+r)+w(1)-a(2))^(-theta-1));
    end
    dist=abs(a2-a);
    if (dist<tol), 
        break
    else
        a2=a;
    end
end
t2=toc;a2
fprintf('the time for Gauss-Siedel algorithm is %8.6f second',t2)
%consumption setting;
disp(' ')
c2=c;
for i=1:1:T-1;
 c2(i)=a2(i)*(1+r)+w(i)-a2(i+1);
 end
 c2(T)=a2(T)*(1+r)+w(T);
 c2
figure
plot(t,c2,'r',t,w,'-*',t,a2,'--o')
title('PS3,Q2,Consumption Saving Problem Linear Gauss-Siedel algorithm')
xlabel('time period')
xlim([1 T])
ylabel('Amount')
legend('c','y','a','Location','NorthEast')


