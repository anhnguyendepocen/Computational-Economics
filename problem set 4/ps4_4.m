clc; clear all; close all;

%% Parameters
bet = 0.96; alp = 0.36; del = 0.06;

%% The steady state
rs = 1/bet - 1 + del;
Ks = (rs/alp)^(1/(alp-1));
ws = (1-alp)*Ks^alp;
Cs = ws + (rs - del)*Ks;
fprintf('Steady state: r = %4.3f, K = %4.3f, w = %4.3f, C = %4.3f. \n', rs, Ks, ws, Cs);

%% Chebyshev nodes and matrix
m = input('m = '); % # nodes (indexed by j)
n = input('n = '); % # Chebyshev functions - 1 (indexed by i)
a = 0.5*Ks; b = 1.5*Ks;

z = -cos( (2*[1:m]-1)/(2*m) * pi )';
K = (z*(b - a) + a + b)/2;

for j = 1:m
    for i = 1:n+1
        T(j,i) = cos((i - 1)*acos(z(j)));
    end
end

%% Chebyshev approximation (when m = n+1!)
tht0 = T\(K.^alp - del.*K); % initial guess: K = Kp
tht = sym('tht', [n+1, 1]);

Kp = K.^alp + (1-del).*K - T*tht;
for j = 1:m
    for i = 1:n+1
        Tp(j,i) = cos((i - 1)*acos( 2*(Kp(j) - a)/(b - a) - 1 ));
    end
end
f = Tp*tht - bet*(1 + alp*Kp.^(alp - 1) - del).*(T*tht);

merr = matlabFunction(f,'Vars',{[tht]});
thts = fsolve(merr,tht0);

syms Ksym; Cpoli = 0;
for i = 1:n+1
    Cpoli = Cpoli + thts(i)*cos((i-1)*acos( 2*(Ksym-a)/(b-a) - 1));
end
Kppoli = Ksym^alp + (1-del)*Ksym - Cpoli;
mCpoli = matlabFunction(Cpoli); mKppoli = matlabFunction(Kppoli);

%% Plot C against K
figure
hold on
set(gcf,'color','w');
xlabel('K_t'); ylabel('C_t');
fplot(mCpoli, [a,b]);
print('CK','-dpng');

%% Plot Kp against K
t = 100; Ksim(1) = 0.5*Ks;
for i = 1:t
    Ksim(i+1) = mKppoli(Ksim(i));
end

figure
hold on
set(gcf,'color','w');
xlabel('K_t'); ylabel('K_{t+1}');
scatter(Ksim([1:end-1]),Ksim([2:end]),'b','fill'); plot(Ksim(1)-0.5:Ksim(end)+0.5,Ksim(1)-0.5:Ksim(end)+0.5);
legend('K_{t+1} = g(K_t)', 'K_{t+1} = K_t', 'Location', 'northwest');
text(Ksim(1)-0.5,Ksim(2)+0.3,'1st period')
text(Ksim(5)-0.5,Ksim(6)+0.3,'5th period')
text(Ksim(10)-0.5,Ksim(11)+0.3,'10th period')
print('KpK','-dpng');

%% Euler equation error
ndraw = 1000;
Krand = 0.5*Ks + Ks.*rand(1,ndraw);

err = zeros(1,ndraw);
for i = 1:ndraw
    err(i) = ( mCpoli( mKppoli(Krand(i)) ) - bet*(1+alp*mKppoli(Krand(i))^(alp-1)-del)*mCpoli(Krand(i)) )/mCpoli(Krand(i));
end

fprintf('The maximal and average errors are %9.8f, %9.8f. \n', max(abs(err)), mean(abs(err)));