% Finds the energy and full wavefunction for a particle from the point 
% of view of an observer in a rotating frame of reference. The 
% wavefunction is correct to order v/c and the spin is in the same 
% direction as the rotation (probably). This is evaluated at a 
% particular time t. This is for a negative energy particle.
% This uses equations (343) and (356)-(361) in my notes in "Dirac 
% equation in a rotating frame of reference 3".
t=2;
k=1.0;
mu=1.5;
om=0.01;
m=1.0;
c=1;
hbar=1;
mf=(m*c^2+hbar*om/2);
E=mu*hbar*om-sqrt(mu^2*hbar^2*om^2+hbar^2*c^2*k^2+mf^2)
X=1;
ef1=E+m*c^2+hbar*om/2;
ef2=E-m*c^2-hbar*om/2;
ef3=E^2-hbar^2*c^2*k^2-(2*mu+2)*hbar*om*E-mf^2;
W=-2*hbar^2*om*c*k*X/ef3;
B=(E*W-hbar*c*k*X)/ef2;
A=(hbar*c*k*W+2*mu*hbar*om*X)/ef2;
rmin=0;
rmax=10;
dr=0.1;
nr=round(1+(rmax-rmin)/dr);
r=rmin:dr:rmax;
z=k*r;
z1=om*z/(k*c);
for ir=1:nr
    psia1(ir)=A*besselj(mu+1,r(ir))+z1(ir)*B*besselj(mu,z(ir));
    psia4(ir)=W*besselj(mu,z(ir))+z1(ir)*X*besselj(mu+1,z(ir));
end
phimax=4*pi;
phimin=0.0*pi;
dphi=0.1*pi;
nphi=1+(phimax-phimin)/dphi;
phi=phimin:dphi:phimax;
et=exp(-1i*E*t/hbar)
for iphi=1:nphi
    Phi(iphi)=exp(1i*(mu+1/2)*phi(iphi));
end
for ir=1:nr
    phi1(ir)=psia1(ir);
    phi4(ir)=-1i*psia4(ir);
    for iphi=1:nphi
        psi1(iphi,ir)=phi1(ir)*Phi(iphi)*et;
        psi4(iphi,ir)=phi4(ir)*Phi(iphi)*et;
    end
end
for ir=1:nr
    pstry1(ir)=psi1(5,ir);
    pstry4(ir)=psi4(5,ir);
end
contour(r,phi,psi1,100)




    
