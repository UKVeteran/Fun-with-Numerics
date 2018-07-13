simin=12.000;
simax=16.000;
dsi=0.01;
nsi=1+(simax-simin)/dsi;
si=simin:dsi:simax;
srmin=0.0;
srmax=1.0;
dsr=0.05;
nsr=1+(srmax-srmin)/dsr;
sr=srmin:dsr:srmax;
for is1=1:nsr
    for is2=1:nsi
        s=sr(is1)+1i*si(is2);
        E=-1i*(s-0.5);
        t1=-E*log(pi)/2;
        t2=-0.5*1i*(log(gammai(s/2))-conj(log(gammai(s/2))));
        theta=t1+t2;
        e1=exp(1i*theta);
        e2=exp(-1i*theta);
        g(is2,is1)=abs(0.5*(e1+e2));
        %g(is2,is1)=abs(cos(theta));
    end
end
contour(sr,si,g,200)