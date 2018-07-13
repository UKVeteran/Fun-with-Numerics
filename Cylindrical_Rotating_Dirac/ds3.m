% Finds the energy and radial wavefunction for a particle from the point 
% of view of an observer in a rotating frame of reference. The 
% wavefunction is correct to order v/c and the spin is in the same 
% direction as the rotation (probably). This is a negative energy solution.  
% This uses equations (343) and (356)-(361) in my notes in "Dirac 
% equation in a rotating frame of reference 3".
k=1.0;
mu=1.5;
om=0.01;
m=1.0;
c=1;
hbar=1;
mf=(m*c^2+hbar*om/2);
E=mu*hbar*om-sqrt(mu^2*hbar^2*om^2+hbar^2*c^2*k^2+mf^2)
Xn=1;
ef1=E+m*c^2+hbar*om/2;
ef2=E-m*c^2-hbar*om/2;
ef3=E^2-hbar^2*c^2*k^2-(2*mu+2)*hbar*om*E-mf^2;
Wn=-2*hbar^2*om*c*k*Xn/ef3;
Bn=(E*Wn-hbar*c*k*Xn)/ef2;
An=(hbar*c*k*Wn+2*mu*hbar*om*Xn)/ef2;
rmin=0;
rmax=10;
dr=0.1;
nr=round(1+(rmax-rmin)/dr);
r=rmin:dr:rmax;
z=k*r;
z1=om*z/(k*c);
A=An/An
B=Bn/An
C=Cn/An
D=Dn/An
for ir=1:nr
    psi1(ir)=A*besselj(mu+1,r(ir))+z1(ir)*B*besselj(mu,z(ir));
    psi2(ir)=W*besselj(mu,z(ir))+z1(ir)*X*besselj(mu+1,z(ir));
end
plot(z,psi1,z,psi2)



    
