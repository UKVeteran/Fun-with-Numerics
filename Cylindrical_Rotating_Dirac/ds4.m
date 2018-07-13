% Finds the energy and radial wavefunction for a particle from the point 
% of view of an observer in a rotating frame of reference. The 
% wavefunction is correct to order v/c and the spin is in the opposite 
% direction as the rotation (probably). This is a negative energy solution.   
k=1.0;
mu=1.5;
om=0.01;
m=1.0;
c=1;
hbar=1;
mf=(m*c^2+hbar*om/2);
E=(mu+1)*hbar*om-sqrt((mu+1)^2*hbar^2*om^2+hbar^2*c^2*k^2+mf^2)
B=1;
ef1=E+m*c^2+hbar*om/2;
ef2=E-m*c^2-hbar*om/2;
ef3=E^2-hbar^2*c^2*k^2-2*mu*hbar*om*E-mf^2;
A=2*hbar^2*om*c*k*B/ef3;
X=(E*A-hbar*k*c*B)/ef1;
W=(hbar*c*k*A+(2*mu+2)*hbar*om*B)/ef1;
rmin=0;
rmax=10;
dr=0.1;
nr=round(1+(rmax-rmin)/dr);
r=rmin:dr:rmax;
z=k*r;
z1=om*z/(k*c);
for ir=1:nr
    psi1(ir)=A*besselj(mu+1,r(ir))+z1(ir)*B*besselj(mu,z(ir));
    psi2(ir)=W*besselj(mu,z(ir))+z1(ir)*X*besselj(mu+1,z(ir));
end
plot(z,psi1,z,psi2)



    
