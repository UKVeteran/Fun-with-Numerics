function yint=simp(dx,n,x,y)
% code to perform a Simpson's rule integration of the array y. x contains n
% equally spaced values of the variable and y contains the integrand at
% these values of the variable.
sumo=0;
for i=3:2:n-2;
    sumo=sumo+y(i);
end
sume=0;
for i=2:2:n-1
    sume=sume+y(i);
end
yint=dx*(y(1)+y(n)+4*sume+2*sumo)/3.;
return
end

