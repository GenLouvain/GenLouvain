function [S,Q,n_it]=iterated_genlouvain(B,limit,verbose,randord,randmove,S0,postprocessor,recursive)

%set default for maximum size of modularity matrix
if nargin<2
    limit = [];
end

%set level of reported/displayed text output
if nargin<3
    verbose = [];
end


%set randperm- v. index-ordered
if nargin<4
    randord = [];
end


%set move function (maximal (original Louvain) or random improvement)
if nargin<5
    randmove=[];
end

% set initial partition
if nargin<6
    S0=[];
end

% set postprocessing function
if nargin<7
    postprocessor=[];
end

% set recursive option
if nargin<8
    recursive=[];
end

S_old=[];
n_it=1;

[S,Q]=genlouvain(B,limit,verbose,randord,randmove,S0,postprocessor,recursive);

Q_old=-inf;
while (~isequal(S,S_old)) && Q-Q_old > 8*eps
    S_old=S;
    Q_old=Q;
    [S,Q]=genlouvain(B,limit,verbose,randord,randmove,S,postprocessor,recursive);
    fprintf('improvement in modularity: %f\n',Q-Q_old)
    n_it=n_it+1;
end

end
