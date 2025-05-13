function [xsol,time]=implicitEuler_SO_ROM(M,D,B,x0,f,del,u,k,V)
%xsol=zeros(n,length(tSim)); % preallocate
zsol(:,1)=V'*x0;
zsol(:,2)=V'*x0;
Mr=V'*M*V;

aa=V'*del*D;
bb=V'*del^2;
dd=V'*del^2;
cc=2*V';
tic
for i=1:k-2
 %zsol(:,i+2)=(Mr)\(V'*(2*V*zsol(:,i+1)-V*zsol(:,i) - del*D*(V*zsol(:,i+1)-V*zsol(:,i)) + del^2*(f(V*zsol(:,i))+B*u(i)))); %%Euler implicit      
xx=V*zsol(:,i); 
xxx=V*zsol(:,i+1);
zsol(:,i+2)=   cc*xxx - V'*xx - aa*(xxx-xx) + bb*(f(xx)) + dd*B*u(i); %%Euler implicit      

end
time=toc;
xsol=V*zsol;
%ysol=C*xsol./n;
end