% This code is to reproduce figure 1 in Pragya and Shukla. It uses
% equations (76) and (77) from my notes on curl forces.
t=linspace(0,20,200);
nt=length(t);
A=1.9;
D=2.2;
x0=.5;
y0=.5;
vx0=.5;
vy0=.5;
M=12;
N=13;
gam=sqrt((A+D)/(M^2+N^2));
for it=1:nt
    f1=(A*M^2-D*N^2)*cos(gam*M*t(it))+(D*M^2-A*N^2)*cos(gam*N*t(it));
    f2=(A+D)*(M^2-N^2);
    f3=cos(gam*N*t(it))-cos(gam*M*t(it));
    f4=gam^2*(M^2-N^2);
    f5=(A*M^2-D*N^2)*N*sin(gam*M*t(it))+(D*M^2-A*N^2)*M*sin(gam*N*t(it));
    f6=gam*(A+D)*M*N*(M^2-N^2);
    f7=M*sin(gam*N*t(it))-N*sin(gam*M*t(it));
    f8=gam^3*M*N*(M^2-N^2);
    g1=(D*M^2-A*N^2)*(A*M^2-D*N^2)*(cos(gam*N*t(it))-cos(gam*M*t(it)));
    g2=(A+D)*(M^4-N^4);
    g3=(A*M^2-D*N^2)*cos(gam*N*t(it))+(D*M^2-A*N^2)*cos(gam*M*t(it));;
    g4=f2;
    g5=(D*M^2-A*N^2)*(A*M^2-D*N^2)*(M*sin(gam*N*t(it))-N*sin(gam*M*t(it)));
    g6=gam*(A+D)*M*N*(M^4-N^4);
    g7=(A*M^2-D*N^2)*M*sin(gam*N*t(it))+(D*M^2-A*N^2)*N*sin(gam*M*t(it));
    g8=gam*(A+D)*M*N*(M^2-N^2);
    x(it)=x0*f1/f2+y0*f3/f4+vx0*f5/f6+vy0*f7/f8;
    y(it)=x0*g1/g2+y0*g3/g4+vx0*g5/g6+vy0*g7/g8;
end
plot(x,y)
    
    
