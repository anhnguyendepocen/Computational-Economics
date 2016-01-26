clear all; close all; clc;

%% question (1)
data    = xlsread('MRW92QJE-data.xls');
data       =  data(:, [1, 3:end]);

% get rid of countries with missing data
D       = isnan(data);
   Dsum = sum(D,2);

 Dx = logical(Dsum == 0);
    data = data(Dx,:);

%% question (2)
%Generate sub samples
DN = logical(data(:,2) == 1);
data1 = data(DN,:); % non-oil

DI = logical(data(:,3) == 1);
data2 = data(DI,:); % intermediate

DO = logical(data(:,4) == 1);
data3 = data(DO,:); % oecd

%% question (3)
%% Regression: Sub-Sample Non-oil

% regression:   dependent variable: difference in log gdp
%               independent variables: constant, log(y60), log(i/gdp),
%               log(n+g+d), log(school)

[rows1, cols1]  = size(data1);
Y1 = log(data1(:,6)) - log(data1(:,5));
X1 = [log(data1(:,5)), log(data1(:,9)/100), log(data1(:,8)/100+0.05), log(data1(:,10)/100)];
mdl1 = fitlm(X1,Y1)    
    
    
   
%% Regression: Sub-Sample Intermediate

% regression:   dependent variable: difference in log gdp
%               independent variables: constant, log(y60), log(i/gdp),
%               log(n+g+d), log(school)

[rows2, cols2]  = size(data2);
Y2 = log(data2(:,6)) - log(data2(:,5));
X2 = [log(data2(:,5)), log(data2(:,9)/100), log(data2(:,8)/100+0.05), log(data2(:,10)/100)];
    
mdl2 = fitlm(X2,Y2)
 

%% Regression: Sub-Sample OECD

% regression:   dependent variable: difference in log gdp
%               independent variables: constant, log(y60), log(i/gdp),
%               log(n+g+d), log(school)

[rows3, cols3]  = size(data3);
Y3 = log(data3(:,6)) - log(data3(:,5));
X3 = [log(data3(:,5)), log(data3(:,9)/100), log(data3(:,8)/100+0.05), log(data3(:,10)/100)];

mdl3 = fitlm(X3,Y3)    


%% question (4) Implied Speed of Convergence: see MRW, p.423, eq. (16): use result from coefficient on ln(Y60) (MRW say so on p. 428)
X1 = [ones(rows1,1), log(data1(:,5)), log(data1(:,9)/100), log(data1(:,8)/100+0.05), log(data1(:,10)/100)];
X2 = [ones(rows2,1), log(data2(:,5)), log(data2(:,9)/100), log(data2(:,8)/100+0.05), log(data2(:,10)/100)];
X3 = [ones(rows3,1), log(data3(:,5)), log(data3(:,9)/100), log(data3(:,8)/100+0.05), log(data3(:,10)/100)];


beta1= X1\Y1;
beta2= X2\Y2;
beta3= X3\Y3;
lambdaa = -log([beta1(2,1), beta2(2,1), beta3(2,1)]+1)/25