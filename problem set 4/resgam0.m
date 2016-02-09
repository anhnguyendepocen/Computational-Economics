function res=resgam0(gam,prob,r,rf);

res=prob'*((1+r).^(-gam).*(r-rf));
