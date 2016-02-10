function util=f_util(gam,c)

epsi=1.0e-04;

if (abs(gam-1.0)<epsi),
    util=log(c);
else
    temp=1.0-gam;
    util=1.0./temp*(c.^temp-1.0);
end

end
