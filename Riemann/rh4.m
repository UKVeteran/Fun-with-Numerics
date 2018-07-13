clear all
smin=-10.0;
smax=10.0;
t=0.7;
ds=0.2;
ns=1+(smax-smin)/ds;
s=smin:ds:smax;
for is=1:ns
    z(is)=0;
    y(is)=zeta(s(is)+1i*0.3);
end
plot(s,real(y),s,imag(y),s,z)

   