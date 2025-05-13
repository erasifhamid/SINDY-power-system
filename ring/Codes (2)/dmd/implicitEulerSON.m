function [xsol,time]=implicitEulerSON(Mr,Dr,Br,x0,fr,Phi,del,u,k)
xsol(:,1)=x0;
xsol(:,2)=x0;
tic
for i=1:k-2
      xx=Phi'*xsol(:,i); 
      xxx=Phi'*xsol(:,i+1);
xsol(:,i+2)=Mr\(2*xxx-xx-del*Dr*(xxx-xx)+del^2*(fr(xx)+Br*u(i)));  %%Euler implicit      
end
time=toc;
end