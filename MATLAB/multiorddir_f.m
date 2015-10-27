function [B,twom]=multiorddir_f(A,gamma,omega)

if nargin<2||isempty(gamma)
    gamma=1;
end

if nargin<3||isempty(omega)
    omega=1;
end

N=length(A{1});
T=length(A);

for i=1:T
    m(i)=sum(A{i}(:));
    
end
A=blkdiag(A{:});
kout=sum(A,1);
koutmat=sparse(1:(N*T),kron(1:T,ones(1,N)),kout);
kin=sum(A,2);
kinmat=sparse(1:(N*T),kron(1:T,ones(1,N)),kin);
A=(A+A')./2;
A=A+omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);

B=@(i) A(:,i)-gamma.*(kout(i).*kinmat(:,ceil(i./(N+eps)))+kin(i).*koutmat(:,ceil(i./(N+eps))))./(2*m(ceil(i./(N+eps))));

twom=sum(m)+omega*2*N*(T-1);
end
