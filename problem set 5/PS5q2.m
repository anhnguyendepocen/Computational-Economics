clc; close all; clear all;
%Parameters setting
r_f=0.02;ru=0.5;mu=[0.04,0.06]';delta_sd=[0.1,0.2];w0=100;wmin=0;gemma=2;a0=6;
rng(321);%set the seed
cov_m=diag(delta_sd);cov=sqrt(ru)*delta_sd(1)*delta_sd(2);
cov_m(2)=cov;cov_m(3)=cov;
r=mu+cov_m*rand(2,1);
%Approximate the expectation by Gauss-Hermite integration using m = 7
%nodes by the MF-function qnwnorm.
m=7;n=[m,m];
[x,w] = qnwnorm(n,mu',cov_m);
%E=w'*f(x),e.g.w'*exp(x(:,1)+x(:,2))
syms a a1 a2 y z1 z2
y=-w'.*1/(1-gemma)*((1+r_f+a*(x(:,1)-r_f)+(1-a)*(x(:,2)-r_f))*w0-wmin).^(1-gemma);
df = diff(y,a);
mdf = matlabFunction(df);
astar1 = fzero(mdf,1);%constant, or to say robost
f=matlabFunction(y);
astar2=fminunc(f,a0);%sensitive to the initial guess
%alternatively to use the FOC we derive to solve the problem
z1=w'*(((1+r_f+a1.*(x(:,1)-r_f)+(1-a1).*(x(:,2)-r_f))*w0-wmin).^(-gemma).*w0.*(x(:,1)-r_f));
z2=w'*(((1+r_f+(1-a2)*(x(:,1)-r_f)+a2*(x(:,2)-r_f))*w0-wmin).^(-gemma).*w0.*(x(:,2)-r_f));
f1=matlabFunction(z1);
f2=matlabFunction(z2);
a1_s1=fzero(f1,a0);
a2_s1=fzero(f2,a0);
right=logical(a1_s1+a2_s1==1);%the answer is wrong, the same if change fsolve to fzero

%So we apply the first method to do replication
wm=(0:10:50);astar=zeros(length(wm),1);
for i=1:1:length(wm)
y=-w'.*1/(1-gemma)*((1+r_f+a*(x(:,1)-r_f)+(1-a)*(x(:,2)-r_f))*w0-wm(i)).^(1-gemma);
df = diff(y,a);
mf = matlabFunction(y);
mdf = matlabFunction(df);
astar(i) = fzero(mdf,1);    
end
plot(wm,astar)

%Homotopy in Lecture 4
%Actually should solve the question by self designed iteration method first
as=a0;itmax=10000;
while i<=itmax
asi=as+w'*(((1+r_f+as.*(x(:,1)-r_f)+(1-as).*(x(:,2)-r_f))*w0-wmin).^(-gemma).*w0.*(x(:,1)-r_f));
%for the cases of two risky assets satisfying two FOCs can be shorten as
%this condition only
error=abs(asi-as);
    if error<sqrt(eps)
    break
    else
        as=asi;i=i+1;
    end
end
asi1=as;
% the result is also robust for guess between 3 and 13.3 so do the loop as well
astar_it=zeros(length(wm),1);as=a0;
for i=1:1:length(wm)
while j<=itmax*10
asi=as+w'*(((1+r_f+as.*(x(:,1)-r_f)+(1-as).*(x(:,2)-r_f))*w0-wm(i)).^(-gemma).*w0.*(x(:,1)-r_f));
error=abs(asi-as);
    if error<sqrt(eps)
    break
    else
    as=asi;j=j+1;
    end
end  
astar_it(i)=as;j=0;as=a0;
end
%it takes around 30 seconds
figure
plot(wm,astar_it)
%Real Homotopy code at last and use symple rule

as=a0;
t=0;
while t<1
while i<=itmax*10

ast=(1-t)*as+t.*w'*(((1+r_f+as.*(x(:,1)-r_f)+(1-as).*(x(:,2)-r_f))*w0-wmin).^(-gemma).*w0.*(x(:,1)-r_f));
error=abs(ast-as);
if error<sqrt(eps)
break
else
    as=ast;i=i+1;
end

end
t=t+0.01;i=0;

end
% take around 30 seconds and use out all steps, 10*i doesn't matter
% also roburt answer,do replication with different asset, too
asi2=as;
syms A Y
A0=[1,1]';
Y=-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wmin).^(1-gemma);
Df = diff(Y,A);
mDf = matlabFunction(Df);
Astar1 = fsolve(mDf,A0);%no solution
Astar2=fsolve(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wmin).^(1-gemma),A0);



%astar, astar1, astar_it, asi1 were thought to be correct answer, but none
%of them is correct
