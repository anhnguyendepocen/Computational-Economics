close all; 
clear all 
clc;
disp(' ')
disp('PS1 Q2')
disp(' ')
%% load data
% read the xls file first to get a broad picture of data,i.e. time length
%open the data manually by importing data or first to change the current folder where the 
%file is and then do the following step
data=xlsread('data/OECD-Germany_Greece_GDP_Linux.xls');
Y_Ge=data(1,:)';Y_Gr=data(2,:)';LY_Ge=log(Y_Ge);LY_Gr=log(Y_Gr);
t=length(data);timeline=linspace(1995.01,2014.04,t)';% set the time 
%% HP filter
lambda=1600; % which is actully the default setting
LY_Ge_HP_T = hpfilter(LY_Ge,lambda);LY_Gr_HP_T = hpfilter(LY_Gr,lambda);
%% OLS filter
tvar=(0:t-1)'; X=[ones(t,1),tvar];% create variable
beta_Ge=X\LY_Ge;beta_Gr=X\LY_Gr;
disp(['the OLS estimator for Germany is beta0=' num2str(beta_Ge(1,1)) ' and beta1=' num2str(beta_Ge(2,1))])
disp(['the OLS estimator for Greece is beta0=' num2str(beta_Gr(1,1)) ' and beta1=' num2str(beta_Gr(2,1))])
LY_Ge_OLS=X*beta_Ge;LY_Gr_OLS=X*beta_Gr;
%% Output Gap
Y_Ge_HP=exp(LY_Ge_HP_T);Y_Gr_HP=exp(LY_Gr_HP_T);Y_Ge_OLS=exp(LY_Ge_OLS);Y_Gr_OLS=exp(LY_Gr_OLS);
G_Y_Ge_HP=(Y_Ge-Y_Ge_HP)./Y_Ge_HP;G_Y_Gr_HP=(Y_Gr-Y_Gr_HP)./Y_Gr_HP;
G_Y_Ge_OLS=(Y_Ge-Y_Ge_OLS)./Y_Ge_OLS;G_Y_Gr_OLS=(Y_Gr-Y_Gr_HP)./Y_Gr_OLS;
%% Plot
LY=[LY_Ge,LY_Gr];LY_HP=[LY_Ge_HP_T,LY_Gr_HP_T];LY_OLS=[LY_Ge_OLS,LY_Gr_OLS];
G_Y_HP=[G_Y_Ge_HP,G_Y_Gr_HP];G_Y_OLS=[G_Y_Ge_OLS,G_Y_Gr_OLS];
zeroline=zeros(t,1);
for n=1:2
figure% open a new figure
plot(timeline, LY(:,n), timeline, LY_HP(:,n), timeline, LY_OLS(:,n))
xlim([min(timeline),max(timeline)]);
xlabel('time');
ylabel('log y_t');
legend('Real Data', 'HP-Trend', 'Linear Trend', 'Location','SouthEast')
if n==1
    title('log GDP and log trend, Germany')
else
    title('log GDP and log trend, Greece')
end
end
for n=1:2
figure% open a new figure
plot(timeline, G_Y_HP(:,n), timeline, G_Y_OLS(:,n),timeline, zeroline)
xlim([min(timeline),max(timeline)]);
xlabel('time');
ylabel('Output Gap');
legend('HP-Trend', 'Linear Trend', 'Location','SouthWest')
if n==1
    title('Output Gap, Germany')
else
    title('Output Gap, Greece')
end
end
