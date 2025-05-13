% GramSchmidt Algorithm
function [Q,R] = GramSchmidt(A)
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j=1:n
    v=A(:,j);        % v begins as column j of A 
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j); %modify A(:,j) to v for more accuracy
        v=v-R(i,j)*Q(:,i);     % subtract the projection (q_i^T*a_j)qi=(q_i)
    end                      
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end
        