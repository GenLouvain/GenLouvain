function B=categorical_coupling(B,omega)
% convert modularity matrices for each layer (given by cellarray `B` into
% multilayer modularity matrix with categorical coupling of strength
% `omega`

T=length(B);
N=length(B{1});
all2all = N*[(-T+1):-1,1:(T-1)];

B=spblkdiag(B{:})+omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);

end
