% problem 1
function problem1

close all

global fspace alpha eta

% -------------------------------------------------------------------------
% test of MF toolbox, see p. 140:
n = 10;
alpha = 2;
fspace = fundefn('cheb',n,-1,1);
c=funfitf(fspace,@myfunc,alpha);
x=nodeunif(n,-1,1);
y=funeval(c,fspace,x);
plot(x,y-myfunc(x,alpha));
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% a math toy:
type='cheb';
alpha=25;

for n=[5,9,13,17],
    disp(' ');
    fspace = fundefn(type,n,-1,1);

    % first: equidistant nodes:
    x=nodeunif(n,-1,1);
    B=funbas(fspace,x);
    y=feval(@mathplay,x,alpha);
    c=B\y;
    x=nodeunif(1001,-1,1);
    y=funeval(c,fspace,x);

    figure;
    res=y-mathplay(x,alpha);
    plot(x,res,'b-');
    disp(['maximum error for cheb with equidistant nodes: ', num2str(max(res))]);

    % now Chebychev nodes:
    % notice: two equivalent formulations:
    % Alternative I:
    c=funfitf(fspace,@mathplay,alpha);

    % Alternative II:
    x=funnode(fspace);
    B=funbas(fspace,x);
    y=feval(@mathplay,x,alpha);
    c=B\y;

    x=nodeunif(1001,-1,1);
    y=funeval(c,fspace,x);

    res=y-mathplay(x,alpha);
    hold on; plot(x,res, 'r--');
    disp(['maximum error for cheb with cheb nodes: ', num2str(max(res))]);

    % now splines with equidistant nodes:
    fspace = fundefn('spli',n,-1,1);
    c=funfitf(fspace,@mathplay,alpha);
    x=nodeunif(1001,-1,1);
    y=funeval(c,fspace,x);

    res=y-mathplay(x,alpha);
    hold on; plot(x,res, 'g-.');
    disp(['maximum error for spline with equidistant nodes: ', num2str(max(res))]);
    disp(' ');
end;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% test collocation of MF, see p. 143
alpha=1.0;
eta=1.5;
n=25;
nlarge=501;
a=1;
b=10;
fspace=fundefn('cheb',n,a,b);
p=funnode(fspace);

% coefficients:
c=ones(n,1);
c=broydn(@myresid,c,1e-6,0,0,p);

% residuals on large grid:
plarge=nodeunif(nlarge,a,b);
res=feval(@myresid,c,plarge);
figure;
plot(plarge,res);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% here's simple adaptive method

% coefficients/range
alpha=1.0;
eta=1.5;
a=1;
b=2.2;

tol=1e-5;
nmin=5;
nstp=4;
nmax=1001;
for t=1:2,
    if t==1,
        type='cheb';
    else
        type='spli';
    end;
    
    fspace=fundefn(type,nmin,a,b);
    p=funnode(fspace);
    for i=nmin:nstp:nmax,
        c=ones(i,1);
        c=broydn(@myresid,c,1e-6,0,0,p);

        % error evaluation
        if i<nmax,
            p=nodeunif(i,a,b);
            res=feval(@myresid,c,p);
            fspace=fundefn(type,i+nstp,a,b);
            if max(abs(res))<tol,
                n=i+nstp;
                p=funnode(fspace);
                c=ones(n,1);
                c=broydn(@myresid,c,1e-6,0,0,p);
                disp(['maximum nodes required: ', num2str(n)]);

                pnew=nodeunif(n,a,b);
                res=feval(@myresid,c,pnew);
                
                % plots:
                figure;
                plot(p,abs(res));
                q=funeval(c,fspace,p);
                figure;
                plot(p,q);
            
                break;
            else,
                p=funnode(fspace);
            end;
        else,
            error('search for optimal n not successful');
        end;
    end;

    if t==1,
        pold=p;
        qold=q;
    end;
end;
figure;
plot(pold,qold,'r--',p,q,'b-');

end
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
function res=myresid(c,p);

global fspace alpha eta

q=funeval(c,fspace,p);
res=p-(q.*p.^(eta+1)./eta)-alpha*sqrt(q)-q.^2;

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function f=myfunc(x,alpha);

f=exp(-alpha*x);

end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
function f=mathplay(x,alpha);

f=1./(1+alpha*x.^2);

end
% -------------------------------------------------------------------------

