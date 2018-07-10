s0=1/2;
ntot=5000;
for is=1:400
    sum=0;
    s=s0+1i*is/100;
for n=1:ntot;
    term=1/n^s;
    sum=sum+term;
end
n2(is)=is;
z(is)=sum;
end
plot(n2,z)
