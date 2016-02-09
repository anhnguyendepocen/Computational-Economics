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

syms A 
A0=[1.2,1]';%very sensitive to initial guess
Astar=fsolve(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wmin).^(1-gemma),A0);

syms a y
y=-w'.*1/(1-gemma)*((1+r_f+a*(x(:,1)-r_f)+(1-a)*(x(:,2)-r_f))*w0-wmin).^(1-gemma);
df = diff(y,a);
mdf = matlabFunction(df);
astar = fzero(mdf,1);%additional constrain with a(1)+a(2)=1

wm=(0:10:50);astar1=zeros(2,length(wm));A1=A0;
for i=1:1:length(wm)
astar1(:,i) = fsolve(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wm(i)).^(1-gemma),A1); 
A1=astar1(:,i);
end
plot(wm,astar1)

% syms y
% a = sym('a', [2 1]);
% y=-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*a)*w0-wmin).^(1-gemma);
% f=matlabFunction(y);
% astar2=fminunc(f,A0);
options = optimoptions('fminunc');
options = optimoptions(options,'Algorithm', 'quasi-newton');
options = optimoptions(options,'MaxFunEvals', 400);
Astar2=fminunc(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wmin).^(1-gemma),A0,options);

wm=(0:10:50);astar2=zeros(2,length(wm));A1=A0;
for i=1:1:length(wm)
astar2(:,i) = fminunc(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wm(i)).^(1-gemma),A1,options); 
A1=astar2(:,i);
end
figure
plot(wm,astar2)

options = optimoptions('fmincon');
options = optimoptions(options,'Algorithm', 'sqp');
Astar3=fmincon(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wmin).^(1-gemma),A0,[],[],[],[],[0;0],[1;1],[],options);
Astar3c=fmincon(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wmin).^(1-gemma),A0,[],[],[1,1],[1],[0;0],[1;1],[],options);
wm=(0:10:50);astar3=zeros(2,length(wm));A1=A0;
for i=1:1:length(wm)
astar3(:,i)=fmincon(@(A)-w'.*1/(1-gemma)*((1+r_f+(x-r_f)*A)*w0-wm(i)).^(1-gemma),A0,[],[],[],[],[0;0],[1;1],[],options);
A1=astar3(:,i);
end
figure
plot(wm,astar3)

optset('ncpsolve','type','smooth');
% USAGE
%     optset(funcname,optname,optvalue)
%   INPUTS
%     funcname : name of function
%     optname  : name of option
%     optval   : option value
alpha0=ones(length(mu),1);
disp(' ');
disp('constrained problem:');
a1=zeros(2,length(wm));
for wmin=0:10:50,
    i=1;
    alph=ncpsolve(@ncpalphres,zeros(length(mu),1),ones(length(mu),1),alpha0,w0,wmin,r_f,x,w);
    alpha0=alph;
    disp(['pf-shares for minimum wealth level ', num2str(wmin), ' are: ', num2str(alph')]); 
    a1(:,i)=alph;i=i+1;
end;
figure
plot(wm,a1(1,:),wm,a1(2,:))