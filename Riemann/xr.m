simin=0.000;
simax=50.000;
ds=0.1;
ns=1+(simax-simin)/ds;
si=simin:ds:simax;
    for is=1:ns
        s=0.5+1i*si(is);
        E=-1i*(s-0.5);
        t1=-E*log(pi)/2;
        t2=-0.5*1i*(log(gammai(s/2))-conj(log(gammai(s/2))));
        theta=t1+t2;
        e1=exp(1i*theta);
        e2=exp(-1i*theta);
        %g(is)=0.5*(e1+e2);
        g(is)=cos(theta);
        z(is)=zeta(s);
    end
    plot(si,real(g))