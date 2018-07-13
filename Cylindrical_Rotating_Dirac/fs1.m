% Finds the energy for a particle from the point 
% of view of an observer in a rotating frame of reference. The 
% wavefunction is correct to order v^2/c^2. This uses equations (362, 
% (363) and (381) of my notes. This is a positive energy particle. 
k=1.0;
mu=1.5;
om=0.0001;
m=1.0;
c=1.0;
hbar=1.0;
hw=hbar*om;
mf=(m*c^2+hbar*om/2);
mom2=hbar^2*c^2*k^2;
emin=m*c*c;
emax=m*c*c+2.0000;
de=0.0001;
e=emin:de:emax;
ne=round(1+(emax-emin)/de);
mf=(m*c^2+hbar*om/2);
for ie=1:ne
    e1=e(ie)^2-mom2-mf^2;
    d1=((e1-2*mu*hbar*om*e(ie))^2)*(e1-2*hbar*om*e(ie));
    d2=8*hbar^4*om^2*c^2*k^2*e(ie)^2;
    d3=4*mu*hbar^2*om^2*e(ie)^2*(e1-2*mu*hbar*om*e(ie));
    det(ie)=d1+d2-d3;
end
det=det
is=0;
for ie=2:ne
    ratio=det(ie)/det(ie-1);
    if ratio < 0
        is=is+1;
        E(is)=e(ie)
    end
end
for ie=2:ne-1
    if det(ie) < det(ie-1)
        if det(ie) < det(ie+1)
            is=is+1
            E(is)=e(ie)
        end
    end
end
