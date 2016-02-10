
function dist=f_dist(utilfun,gam,y1,R,w)

f1=feval(utilfun,gam,y1);
f2=feval(utilfun,gam,R);
f2=w'*f2;

dist=f1-f2;

end