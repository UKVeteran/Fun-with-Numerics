simin=12.000;
simax=16.000;
dsi=0.01;
nsi=round(1+(simax-simin)/dsi);
si=simin:dsi:simax;
srmin=0.000;
srmax=2.000;
dsr=0.05;
nsr=round(1+(srmax-srmin)/dsr);
sr=srmin:dsr:srmax;
for is1=1:nsr
    for is2=1:nsi
        s=sr(is1)+1i*si(is2);
        s2=s/2;
        g1=gammai(s2);
        f1=pi^(s2)/g1;
        sstar=sr(is1)-1i*si(is2);
        sstar2=sstar/2;
        g2=gammai(sstar2);
        f2=pi^((1-s)/2)/g2;
        the(is2,is1)=abs(f1+f2);
    end
end
contour(sr,si,the,200)