function result=value(k,sr,a,b,s,t,n)
result=(((ones(1,n)-sr).*(1-a).*(k.^a)).^(1-s)+b.*(((k.^a).*sr).^(a.*(1-s)).*t))/(1-s)+logical(sr>1)*(-10^7);%1*n dimension
end

       