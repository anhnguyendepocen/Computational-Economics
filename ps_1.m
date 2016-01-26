clc
close all
clear all
format long

%% Question (5)
disp('---------------------------------------------')
disp('Problem Set 1 - Exercise 1 - question (5)')

%% Show convergance Graphically for the convergence case
A=[1 0.5; 1 -1];
y=[3 1]';
x_ini=[0.1 0.1]'; %initial guess
Q=tril(A); 
dx=ones(2,1);
k=0;
q=0:5;
hold on
a=3;b=.5;c=1;d=1;
p = a - b * q;
plot(q,p);
p = c + d * q;
plot(q,p);
title('The case of convergence')
xlabel('q')
ylabel('p')
x=x_ini;
while (dx > sqrt(eps))
k=k+1;
x_new= inv(Q)*y + (eye(size(A)) - inv(Q)*A )*x; % Gauss-Seidel Method
dx=abs(x_new-x);
x=x_new;
plot(x(2),x(1),'bo');
text(x(2),x(1),num2str(k));
if k>1000 
disp('Not converged');
break
end
end
answer=strcat('Equlibrium Price is P*= ',num2str(x_new(1)),' & Equlibrium Quantity is q*= ',num2str(x_new(2)));
disp(answer)

figure;
hold on
A=[1 -1; 1 0.5];
y=[1 3]';
x_ini=[0.1 0.1]'; %initial guess
Q=tril(A);
dx=ones(2,1);
q=-100:1:100;
a=3;b=.5;c=1;d=1;
p = c + d * q;
plot(q,p);
p = a - b * q;
plot(q,p)
title('The case of non-convergence')
xlabel('q')
ylabel('p')
x=x_ini;
k=0;
while (dx > sqrt(eps))
k=k+1;
x_new= inv(Q)*y + (eye(size(A)) - inv(Q)*A )*x; % Gauss-Seidel Method
dx=abs(x_new-x);
x=x_new;
plot(x(2),x(1),'yo');
text(x(2),x(1),num2str(k));
if k>5
disp('Not converged');
break
end
end
hold off




%% Question (6)
disp('---------------------------------------------')
disp('Problem Set 1 - Exercise 1 - question (6)')

x=x_ini;
k=0;
for lambda=0.1:0.1:0.9
while (dx > sqrt(eps))
k=k+1;
x_new= inv(Q)*y + (eye(size(A)) - inv(Q)*A )*x; %Gauss©\Seidel Method
dx=abs(x_new-x);
x=x_new*lambda + x*(1-lambda);
if k>1000
disp('Not converged');
break
end
end
time(round(lambda/.1))=k;
end
[~,I] = min(time);
answer=strcat('Minimum convergence time is for lambda=',num2str(I*.1));
disp(answer)
answer=strcat('Equlibrium Price is P*= ',num2str(x_new(1)),' & Equlibrium Quantity is q*= ',num2str(x_new(2)));
disp(answer)





















