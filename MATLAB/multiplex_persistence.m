function pers=multiplex_persistence(S)

if ~iscell(S)
    S={S};
end

[N,T]=size(S{1});
all2all = N*[(-T+1):-1,1:(T-1)];
A=spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
pers=zeros(length(S),1);
for i=1:length(S)
    
    G=sparse(1:length(S{i}(:)),S{i}(:),1);
    pers(i)=trace(G'*A*G)/(N*T*(T-1));
end

end
