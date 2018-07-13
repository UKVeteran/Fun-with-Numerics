clear all
trmin=0.0000;
trmax=10.000;
timin=0.0;
timax=1.0;
dtr=0.05;
dti=0.05;
ntr=round(1+(trmax-trmin)/dtr);
nti=round(1+(timax-timin)/dti);
tr=trmin:dtr:trmax;
ti=timin:dti:timax;
for itr=1:ntr
    for iti=1:nti
       t(iti,itr)=tr(itr)+1i*ti(iti);
       term1=-t(iti,itr)*log(pi)/2.;
       arg1=1/4+1i*t(iti,itr)/2;
       term2=imag(log(gammai(arg1)));
       theta(iti,itr)=term1+term2;
       eith(iti,itr)=exp(1i*theta(iti,itr));
       arg2=1/4-1i*conj(t(iti,itr))/2;
       s1=gammai(1/4+1i*t(iti,itr)/2);
       s2=gammai(1/4-1i*conj(t(iti,itr))/2);
       trait(iti,itr)=pi^(-1i*t(iti,itr)/2)*sqrt(s1/s2);
       cos2(iti,itr)=(eith(iti,itr)+1/eith(iti,itr))/2;
       phi(iti,itr)=sqrt((1+cos2(iti,itr))/2);
    end
end
for itr=1:ntr
    tt(itr)=trait(21,itr);
    ee(itr)=eith(21,itr);
end
plot(tr,real(ee),tr,real(tt))
%contour(tr,ti,phi,200)
