% Finds the energy and full wavefunction for a particle from the point 
% of view of an observer in a rotating frame of reference. The 
% wavefunction is correct to order v^2/c^2. This uses equations (362), 
% (363) and (373) to (381) of my notes. This is a positive energy particle. 
t=0;
k=1.0;
mu=2.5;
om=0.01;
m=1.0;
c=5.0;
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
ns=0;
for ie=2:ne
    ratio=det(ie)/det(ie-1);
    if ratio < 0
        ns=ns+1;
        E(ns)=e(ie)
    end
end
for ie=2:ne-1
    if det(ie) < det(ie-1)
        if det(ie) < det(ie+1)
            ns=ns+1
            E(ns)=e(ie)
        end
    end
end
for is=1
    et=exp(-1i*E(is)*t/hbar);
    ef1=(E(is)+m*c^2+hbar*om/2);
    ef3=E(is)^2-mom2-2*mu*hbar*om*E(is)-mf^2;
    A=1
    C=E(is)^2*A/(E(is)^2-mom2-2*hbar*om*E(is)-mf^2)
    B=-4*hbar^2*om^2*c*k*C/ef3
    W=(hbar*c*k*A+(2*mu+2)*hbar*om*B)/ef1
    X=(E(is)*A-hbar*c*k*B+2*hbar*om*C)/ef1
    Y=(E(is)*B+hbar*c*k*C)/ef1
end
rmin=0;
rmax=30;
dr=0.1;
nr=round(1+(rmax-rmin)/dr);
r=rmin:dr:rmax;
z=k*r;
z1=om*z/(k*c);
z2=z1.*z1;
for ir=1:nr
    psia1(ir)=A*besselj(mu+1,z(ir))+B*z1(ir)*besselj(mu,z(ir))+C*z2(ir)*besselj(mu+1,z(ir));
    psia4(ir)=W*besselj(mu,z(ir))+X*z1(ir)*besselj(mu+1,z(ir))+Y*z2(ir)*besselj(mu,z(ir));
end
phimax=4*pi;
phimin=0.0*pi;
dphi=0.01*pi;
nphi=round(1+(phimax-phimin)/dphi);
phi=phimin:dphi:phimax;

psi1=zeros(nphi,nr);
psi4=zeros(nphi,nr);
for iphi=1:nphi
    Phi(iphi)=exp(1i*(mu+1/2)*phi(iphi));
end
for ir=1:nr
    phi1(ir)=psia1(ir)*et;
    phi4(ir)=-1i*psia4(ir)*et;
    for iphi=1:nphi
        psi1(iphi,ir)=phi1(ir)*Phi(iphi);
        psi4(iphi,ir)=phi4(ir)*Phi(iphi);
    end
end
for iphi=1:nphi
    pstry1(iphi)=psi1(iphi,20);
    pstry4(iphi)=psi4(iphi,20);
end
plot(phi,pstry1,phi,pstry4)

    
    
    