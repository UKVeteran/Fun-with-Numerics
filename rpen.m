% Calculates something
h=4;
q0=20.0;
N0=4000;
q=q0*sqrt(h);
N=N0*h;
p=linspace(0.0,q,N);
N=length(p);
niter=1000;
A=zeros(length(p),length(p));
V=zeros(N,1);
f=zeros(N,1);
mod=zeros(1,niter);
for ip=1:N
    V(ip)=rand;
end
mod(1)=norm(V);
m=1.0;
c=1.0;
hbar=1.0;
T=25;
eps=sqrt(4*hbar/(m*c^2*T))
eps2=eps^2.0;
r=p./(m*c*eps);
q=p./(m*c);
eta=sqrt(1+eps^2*r.^2);
for ip=1:N
    p1(ip)=sqrt(eta(ip)*(eta(ip)+1.0));
    etap1(ip)=eta(ip)+1.0;
end
for ip=1:N
    for iq=1:ip-1
        nfac=(r(ip)*etap1(iq)+r(iq)*etap1(ip))/(p1(ip)*p1(iq));
        aa=2.0*(eta(ip)-eta(iq))/eps2;
        A(ip,iq)=-nfac*sin(aa)/(pi*aa);
        A(iq,ip)=A(ip,iq);
    end
end
for ip=1:N
    nfacn=2.0*r(ip)*(eta(ip)+1.0);
    nfacd=eta(ip)*(eta(ip)+1.0);
    nfac=nfacn/nfacd;
    A(ip,ip)=-sin(1)*nfac/pi;
end
A=A*r(N)/N;
for ip=1:N
    A(ip,ip)=A(ip,ip)+1.0;
end
%[V,D]=eig(A);
%for ik=1:N
%    zz(ik)=D(ik,ik);
%end
%zz=zz
%stop
for it=2:niter;
    it=it
    Vold=V;
    V=A*V/mod(it-1);
    mod(it)=norm(V);
    V=V;
end
a=0.0;
for ik=1:N
    a=a+V(ik)*Vold(ik)/mod(it-1);
end
a=a-1.0
fileID = fopen('bf.txt','w');
%fprintf(fileID,'%6s %12s\n','p','V(p)');
%for ie=1:N
fprintf(fileID,'%12.8f %12.8f\n',p,V);
%end
fclose(fileID);




plot(q,V)
    
        
    

    
