function [xsol,time]=ROM_DEIM(M,D,B,x0,f,del,u,k,V,int,pre,P)
%xsol=zeros(n,length(tSim)); % preallocate
xsol(:,1)=x0;
xsol(:,2)=x0;

aa=V'*del*D;
bb=del^2*pre;
dd=V'*del^2;
cc=2*V';
tic
for i=1:k-2    
xx=V*xsol(:,i); 
xxx=V*xsol(:,i+1);
xsol(:,i+2)=   cc*xxx - V'*xx - aa*(xxx-xx) + bb*(f(P'*xx)) + dd*B*u(i); %%Euler implicit      
end
time=toc;
%ysol=C*xsol./n;
end