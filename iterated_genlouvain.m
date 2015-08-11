function [S,Q,n_it]=iterated_genlouvain(B,limit,verbose,randord,randmove,S0,postprocessor,recursive)

%set default for maximum size of modularity matrix
if nargin<2||isempty(limit)
    limit = 10000;
end

%set level of reported/displayed text output
if nargin<3||isempty(verbose)
    verbose = 1;
end
if verbose
    mydisp = @(s) disp(s);
else
    mydisp = @(s) disp('');
end

%set randperm- v. index-ordered
if nargin<4||isempty(randord)
    randord = 1;
end
if randord
    myord = @(n) randperm(n);
else
    myord = @(n) 1:n;
end

%set move function (maximal (original Louvain) or random improvement)
if nargin<5||isempty(randmove)
    randmove=false;
end
if randmove
    movefunction='moverand';
else
    movefunction='move';
end

% set initial partition
if nargin<6||isempty(S0)
    S0=[];
end

% set postprocessing function
if nargin<7||isempty(postprocessor)
    % default: do not do anything
    postprocessor=@(S) S;
end

% set recursive option
if nargin<8||isempty(recursive)
    recursive=false;
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
