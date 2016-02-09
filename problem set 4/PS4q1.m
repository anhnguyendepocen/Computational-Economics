clc; close all; clear all;
%% Problem Set 4, Question 1
% USE fspace = fundefn(bastype,n,a,b,order) in M&F
% Setup Function
syms x y
y=1/(1+25*x^2);
f=matlabFunction(y);
% Setup Parameters
bastype='cheb';n=100;%nth-degree Chebychev approximants;
a=-1;b=1;node=(5:2:15);l=length(node);
fspace = fundefn(bastype,n,a,b);%fspace is a structured MATLAB variable containing
%numerous fields of information necessary for forming approximations in the chosen
%function space
c = funfitf(fspace,f);%compute the coefficient vector for the approximant
%that interpolates the function at the standard Chebychev nodes
x_v = nodeunif(1001,-1,1);%simulate real function
y_ch = funeval(c,fspace,x_v);%simlated function
plot(x_v,y_ch-f(x_v));%compare difference
figure;plot(x_v,y_ch)
x_ch = funnode(fspace);
figure;o=zeros(length(x_ch),1);scatter(x_ch,o)

%% To simulate on fix nodes
x_eq=(-1:0.01:1)';%define a vector of equidistant nodes x
B=funbas(fspace,x_eq);%returns the matrix containing the values of the basis functions evaluated at the points x.
y_eq=B*c;%you next calculate the function values at x
c_B=B\y_eq;%get the polynomial coefficients,think as OLS!!
figure;plot(x_eq,y_eq);

figure;
hold on
for i=1:1:l
    x_eq_i=(-1:2/(-1+node(i)):1)';
    B_i=funbas(fspace,x_eq_i);
    y_eq_i=B_i*c;
    plot(x_eq_i,y_eq_i);
end
hold off
%Repeat the exercise using Chebychev nodes.
figure;hold on
for i=1:1:l
fspace_i = fundefn(bastype,node(i),a,b);
x_ch_i = funnode(fspace_i);
c_i = funfitf(fspace_i,f);
y_ch_i = funeval(c,fspace,x_v);
plot(x_v,y_ch_i);
end
hold off

figure;hold on
for i=1:1:l
fspace_i = fundefn(bastype,node(i),a,b);
x_ch_i = funnode(fspace_i);
c_i = funfitf(fspace_i,f);
y_ch_i = funeval(c,fspace,x_ch_i);
plot(x_ch_i,y_ch_i);
end
hold off

figure
for i=1:1:l
fspace_i = fundefn(bastype,node(i),a,b);
x_ch_i = funnode(fspace_i);
o_i=zeros(length(x_ch_i),1);
subplot(l,1,i);
scatter(x_ch_i,o_i)
end

x_vand=vander(x_v);
figure;
for i=1:1:l
    fspace_i = fundefn(bastype,node(i),a,b);
    c_i = funfitf(fspace_i,f);
    x_eq_i=(-1:2/(-1+node(i)):1)';
    B_i=funbas(fspace_i,x_eq_i);
    y_eq_i=B_i*c_i;
    c_B_i=B_i\y_eq_i;
    y_eqv=x_vand(1:node(i),:)'*c_B_i;
    subplot(2,3,i);
    scatter(x_v,y_eqv);
end


