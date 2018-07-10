    % This code is to reproduce figure 1 in Pragya and Shukla. It uses
% equations (76) and (77) from my notes on curl forces.
t=linspace(0,2000,5000);
nt=length(t);
A=2
G=4;
x0=11;
y0=12;
z0=13;
vx0=12.3;    
vy0=32;
vz0=12.12;
M=30;
N=80;
L=40;
m2=-M^2*(A+G)/(M^2+N^2);
m3=-N^2*(A+G)/(M^2+N^2);
d=N^2*m2/L^2;
D=-d;
gam=sqrt((A+G)/(M^2+N^2));
sd=sqrt(D);
for it=1:nt
    f1=(A*M^2-G*N^2)*cos(gam*M*t(it))+(G*M^2-A*N^2)*cos(gam*N*t(it));
    f2=(A+G)*(M^2-N^2);
    f3=cos(gam*N*t(it))-cos(gam*M*t(it));
    f4=gam^2*(M^2-N^2);
    f5=(A*M^2-G*N^2)*N*sin(gam*M*t(it))+(G*M^2-A*N^2)*M*sin(gam*N*t(it));
    f6=gam*(A+G)*M*N*(M^2-N^2);
    f7=M*sin(gam*N*t(it))-N*sin(gam*M*t(it));
    f8=gam^3*M*N*(M^2-N^2);
    g1=(G*M^2-A*N^2)*(A*M^2-G*N^2)*(cos(gam*N*t(it))-cos(gam*M*t(it)));
    g2=(A+G)*(M^4-N^4);
    g3=(A*M^2-G*N^2)*cos(gam*N*t(it))+(G*M^2-A*N^2)*cos(gam*M*t(it));
    g4=f2;
    g5=(G*M^2-A*N^2)*(A*M^2-G*N^2)*(M*sin(gam*N*t(it))-N*sin(gam*M*t(it)));
    g6=gam*(A+G)*M*N*(M^4-N^4);
    g7=(A*M^2-G*N^2)*M*sin(gam*N*t(it))+(G*M^2-A*N^2)*N*sin(gam*M*t(it));
    g8=gam*(A+G)*M*N*(M^2-N^2);
    x(it)=x0*f1/f2+z0*f3/f4+vx0*f5/f6+vz0*f7/f8;
    y(it)=y0*cos(sd*t(it))+vy0*sin(sd*t(it))/sd;
    z(it)=x0*g1/g2+z0*g3/g4+vx0*g5/g6+vz0*g7/g8;
end
plot3(x,y,z)
    
    
