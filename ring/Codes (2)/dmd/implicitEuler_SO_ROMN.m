function [xsol,time]=implicitEuler_SO_ROMN(M,D,B,x0,f,del,u,k,Vp1)
%xsol=zeros(n,length(tSim)); % preallocate;
zsol(:,1)=Vp1'*x0;
zsol(:,2)=Vp1'*x0;
aa=Vp1'*del*D;
bb=Vp1'*del^2;
dd=Vp1'*del^2;
cc=2*Vp1';
tic
for i=1:k-2
 %zsol(:,i+2)=(Mr)\(V'*(2*V*zsol(:,i+1)-V*zsol(:,i) - del*D*(V*zsol(:,i+1)-V*zsol(:,i)) + del^2*(f(V*zsol(:,i))+B*u(i)))); %%Euler implicit      
xx=V*zsol(:,i); 
xxx=V*zsol(:,i+1);
zsol(:,i+2)=   cc*xxx - Vp1'*xx - aa*(xxx-xx) + bb*(f(xx)) + dd*B*u(i); %%Euler implicit      

end
time=toc;
xsol=V*zsol;
%ysol=C*xsol./n;
end