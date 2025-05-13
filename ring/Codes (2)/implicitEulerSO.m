function [xsol,time]=implicitEulerSO(M,D,B,x0,f,del,u,k)
%xsol=zeros(n,length(tSim)); % preallocate
xsol(:,1)=x0;
xsol(:,2)=x0;
tic
for i=1:k-2
 xsol(:,i+2)=M\(2*xsol(:,i+1)-xsol(:,i)-del*D*(xsol(:,i+1)-xsol(:,i))+del^2*(f(xsol(:,i))+B*u(i))); %%Euler implicit      
end
time=toc;
%ysol=C*xsol./n;
end