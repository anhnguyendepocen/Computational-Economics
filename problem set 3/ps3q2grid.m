clc;clear all; close all;
%% Grid Method backward induction
%  do backwards induction starting from period T-1 to obtain the decision 
%  rules in all periods given all possible asset levels
%  suppose you are given all possible a(t), what is the optimal a(t+1)
%% Parameters
beta=0.95;r=0.02;theta=1;T=50;%R=40; if there is a retirement period

%income setting,can be some progress
y=zeros(1,T);y(1)=10;t=(1:1:T);

%form asset grid, size sensitive to the asset and income
%time index for asset is critical, a(t) is the amount you get in the beginning
%of a period, so a(0)=0,a(T+1)=0
a_grid=(-25:0.1:25);%1*N
N=length(a_grid);

%initialize the value function & utility function
v=zeros(T,N);
c=y(T)*ones(1,N)+a_grid;%1*N+1*N=1*N for the last period,atcually doesn't matter here
if theta==1,
    u=log(c);
else
    u=(c.^(1-theta)-1)/(1-theta);
end
ind_c=find(c<=0);
u(ind_c)=-1e10;%penalty to negative consumption
v(T,:)=u(1,:);%T*N,value in period

a_opt(T,:)=zeros(1,N);%T*N,to set up a optimal asset, little tricky in matlab 
%language to create such a matrix, same to code as a_opt=zeros(T,N);
ind_opt(T,:)=ones(1,N);%T*N,to find out in each period, which one is the best choice

for i=(T-1):-1:1 % start to find optimal in T-1 first
    
    % Extremly important!! Imagine a(t+1) as y, a(t) as x
    % Note: The asset level today is in the columns and the choices for 
    % tomorrow are on the rows (therefore transpose the potential choices a_t+1)
    c=y(i)*ones(N,N)+ones(N,1)*a_grid*(1+r)-a_grid'*ones(1,N);
    %c(t)=w(t)+a(t)*(1+r)-a(t+1),n*n+n*n
    if theta==1,
        u=log(c);
    else
        u=(c.^(1-theta)-1)/(1-theta);
    end
    ind_c=find(c<=0);
    u(ind_c)=-1e10;%penalty to negative consumption


    % Note: Columns: state of today, rows: choices for tomorrow (therefore 
    % transpose v_t+1)
    [v(i,:),ind_max]=max(u+beta*v(i+1,:)'*ones(1,N),[],1);
    % returns the largest elements along the dimension of A specified
    %by scalar dim. For example, max(A,[],1) produces the maximum values 
    %along the first dimension of A, column by column.
    %then give this row value to the i period of value funtion
    a_opt(i,:)=a_grid(1,ind_max);
    %determine the optimal asset for period i,
    ind_opt(i,:)=ind_max;
    %as ind_max in the loop, to keep record of very t period

end

%  so far we have finished grid searching

a_sim=zeros(T,1);
c_sim=zeros(T,1);
a_sim(1)=0;% to simulate the savings and consumption paths of an agent with 
%assets equal to zero at the beginning of life 

for i=1:(T-1)
    
    ind_asset=find(a_grid>=a_sim(i,1),1,'first');
    %returns at most the first k indices corresponding to the nonzero
    %entries of X. k must be a positive integer.
    %to find the same asset and locate on the grid
    a_sim(i+1,1)=a_opt(i,ind_asset);
    %as we know the value of a(t) on the grid, we find on the a_opt the
    %optimal value of a(t+1) and give it to a_sim(i+1,1)
    c_sim(i,1)=y(i)+a_sim(i,1)*(1+r)-a_sim(i+1,1);
    %c(t)=w(t)+a(t)*(1+r)-a(t+1),
    
end
c_sim(T,1)=y(T)+a_sim(T,1);%last period

%  plot the consumption and savings paths

plot(t,c_sim,'r',t,y,'-*',t,a_sim,'--o')
title('PS3,Q2,Consumption Saving Problem')
xlabel('time period')
xlim([1 T])
ylabel('Amount')
legend('c','y','a','Location','NorthEast')