close all; 
clear all 
clc;
disp(' ')
disp('PS2 Q4')
disp(' ')
%% 1 
disp('Q1');
sk=0.200;sh=0.200;n=0.010;g=0.015;deltak=0.10;deltah=0.06;z=1;alphak=0.33;alphah=0.33;
params=[sk;sh;n;g;deltak;deltah;z;alphak;alphah;];
k_star=((sk/(n+g+n*g+deltak))^(1-alphah)*(sh/(n+g+n*g+deltah))^alphah)^(1/(1-alphak-alphah));
h_star=((sk/(n+g+n*g+deltak))^alphak    *(sh/(n+g+n*g+deltah))^(1-alphak))^(1/(1-alphak-alphah));
x_star=[k_star; h_star];
disp(['k_star is ' num2str(x_star(1,1)) ' and h_star is ' num2str(x_star(2,1)) ' for given parameters'])
%% 2
f = @(x)ps242(x,params);% return to value and Jacobian of steady state condition given any vector,not a to find a solution. 
%% 3 
disp('Q3 Newton')
tole=eps;told=eps;miter=1000000;I=[tole,told,miter];%set tolerence and maxnumber of iteration
x0=[3,3]';%set initial value of k&h
[x01,fx01,ef01,iter01] = ps243(f,x0,I);
if (~isreal(x01)) % depending on the initial guess, the algorithm may enter negative consumption or capital and return complex values, e.g. x0=[1;1]
        error('solution contains imaginary part, try different starting value');
end
disp(' ')
fprintf('exitflag: %.2d solution: k = % 6.4f, h = %6.4f \n',ef01,x01);
fprintf('number of iterations:%2.0f',iter01);
disp(['  infinity norm: ' num2str(norm(x01-x_star,Inf))])
fprintf('tole:%e told:%e miter:%2.0f',I(1,1),I(1,2),I(1,3));
I_other=[sqrt(I);sqrt(sqrt(I))];
size_I=size(I_other);
for n=1:size_I(1,1)
    [x01,fx01,ef01,iter01] = ps243(f,x0,I_other(n,:));
if (~isreal(x01)) % depending on the initial guess, the algorithm may enter negative consumption or capital and return complex values, e.g. x0=[1;1]
        error('solution contains imaginary part, try different starting value');
end
disp(' ');disp(' ');
fprintf('exitflag: %.2d solution: k = % 6.4f, h = %6.4f \n',ef01,x01);
fprintf('number of iterations:%2.0f',iter01)
disp(['  infinity norm: ' num2str(norm(x01-x_star,Inf))])
fprintf('tole:%e told:%e miter:%2.0f',I_other(n,1),I_other(n,2),I_other(n,3));
disp(' ');
end
%% 4
disp('Q 4 broyden')
[x02,fx02,ef02,iter02] = ps244(f,x0,I);
if (~isreal(x02)) % depending on the initial guess, the algorithm may enter negative consumption or capital and return complex values, e.g. x0=[1;1]
        error('solution contains imaginary part, try different starting value');
end
disp(' ')
fprintf('exitflag: %.2d solution: k = % 6.4f, h = %6.4f \n',ef02,x02);
fprintf('number of iterations:%2.0f',iter02);
disp(['  infinity norm: ' num2str(norm(x02-x_star,Inf))])
fprintf('tole:%e told:%e miter:%2.0f',I(1,1),I(1,2),I(1,3));
disp(' ')
%% 5
disp('Q 5 inverse broyden');
[x03,fx03,ef03,iter03] = ps245(f,x0,I);
if (~isreal(x03)) % depending on the initial guess, the algorithm may enter negative consumption or capital and return complex values, e.g. x0=[1;1]
        error('solution contains imaginary part, try different starting value');
end
disp(' ')
fprintf('exitflag: %.2d solution: k = % 6.4f, h = %6.4f \n',ef03,x03);
fprintf('number of iterations:%2.0f',iter03);
disp(['  infinity norm: ' num2str(norm(x03-x_star,Inf))])
fprintf('tole:%e told:%e miter:%2.0f',I(1,1),I(1,2),I(1,3));
disp(' ')
%% 6 
disp('Q 6 fixed point');
TS=zeros(2,I(1,3));TS(:,1) = x0; 
for j = 1:I(1,3)
    [xp] = ps246(x0,params);
    if norm(x0-xp) <= I(1,1)*(1+norm(xp))
        break;
    else
        TS(:,j+1) = xp; x0 = xp;% save all results
    end
end
iter04 = j; T=TS(:,1:j); clear TS;

figure;
plot((1:iter04),T(1,:),'--k', (1:iter04),T(2,:), '-b')
hold on
line(get(gca,'xlim'),[k_star k_star]); % Adapts to x limits of current axes
line(get(gca,'xlim'), [h_star h_star]); % Adapts to x limits of current axes
hold off
xlabel('Iteration','FontSize',16)
ylabel('paths of k and h','FontSize',16)
legend('k','h', 'Location','Best')
set(gca,'FontSize',16)