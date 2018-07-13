% Finds the energy and radial wavefunction for a particle from the point 
% of view of an observer in a rotating frame of reference. The 
% wavefunction is correct to order v^2/c^2. This uses equations (362, 
% (363) and (382) to (390) of my notes. This is a positive energy particle. 
% This is for a negative energy particle.
k=3.0;
mu=4;
om=0.0000000000001;
m=1.0;
c=1.0;
hbar=1.0;
hw=hbar*om;
mf=(m*c^2+hbar*om/2);
mom2=hbar^2*c^2*k^2;
emin=-m*c*c-2.0000;
emax=-m*c*c;
de=0.0001;
e=emin:de:emax;
ne=round(1+(emax-emin)/de);
mf=(m*c^2+hbar*om/2);
for ie=1:ne
    e1=e(ie)^2-mom2-mf^2;
    d1=(e1-2*(mu+1)*hbar*om*e(ie))^2*(e1+2*hbar*om*e(ie));
    d2=8*hbar^4*om^2*c^2*k^2*e(ie)^2;
    d3=4*(mu+1)*hbar^2*om^2*e(ie)^2*(e1-2*(mu+1)*hbar*om*e(ie));
    det(ie)=d1+d2+d3;
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
for is=1:1
     ef3=E(is)^2-mom2-2*mu*hbar*om*E(is)-mf^2;
    Wn=1;
    Yn=E(is)^2*Wn/(E(is)^2-mom2+2*hbar*om*E(is)-mf^2);
    Xn=(4*hbar^2*om*c*k/(E(is)^2-mom2-2*(mu+1)*hbar*om*E(is)-mf^2))*Yn;
    An=(hbar*c*k*Wn+2*mu*hbar*om*Xn)/ef3;
    Bn=(E(is)*Wn-hbar*c*k*Xn-2*hbar*om*Yn);
    Cn=(E(is)*Xn+hbar*c*k*Yn)/ef3;
end
rmin=0;
rmax=10;
dr=0.1;
nr=round(1+(rmax-rmin)/dr);
r=rmin:dr:rmax;
z=k*r;
z1=om*z/(k*c);
z2=z1.*z1;
A=An/Yn
B=Bn/Yn
C=Cn/Yn
W=Wn/Yn
X=Xn/Yn
Y=Yn/Yn
for ir=1:nr
    psi1(ir)=A*besselj(mu+1,z(ir))+B*z1(ir)*besselj(mu,z(ir))+C*z2(ir)*besselj(mu+1,z(ir));
    psi4(ir)=W*besselj(mu,z(ir))+X*z1(ir)*besselj(mu+1,z(ir))+Y*z2(ir)*besselj(mu,z(ir));
end
plot(z,psi1,z,psi4)
    