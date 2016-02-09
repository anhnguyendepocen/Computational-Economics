function res=alphres(alph,w0,wmin,rf,rnodes,rwghts)

global gemma;

n=length(alph);
res=zeros(n,1);
for i=1:n,
    res(i) = rwghts'*(((w0.*(1+rf+(rnodes-rf)*alph)-wmin).^(-gemma)).*w0.*(rnodes(:,i)-rf));
end;