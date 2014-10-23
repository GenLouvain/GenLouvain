function [B,twom]=multiplex_modularity(A,gamma,omega,wrapper,coupling)

N=length(A{1});
T=length(A);
B=zeros(N*T);
for i=1:T
    ind=(1+(i-1)*N):(i*N);
    [B(ind,ind),twom(i)]=wrapper(A{i},gamma);
end

switch coupling
    case 'ordinal'
        B=B+omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
        twom=sum(twom)+2*N*(T-1)*omega;
        
    case 'categorical'
        B=B+omega*spdiags(ones(N*T,2*T-2),N*[(-T+1):-1,1:(T-1)],N*T,N*T);
        twom=sum(twom)+N*T*(T-1)*omega;
    otherwise
        error('unknown coupling')
end
        


