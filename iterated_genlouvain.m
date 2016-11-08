function [S,Q,n_it]=iterated_genlouvain(B,limit,verbose,randord,randmove,S0,postprocessor)
% Optimise modularity-like quality function by repeatedly applying GenLouvain
%
% Version:
% Date:

%set default for maximum size of modularity matrix
if nargin<2
    limit = [];
end

%set level of reported/displayed text output
if nargin<3||isempty(verbose)
    verbose = false;
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

% verbose output switch
if verbose
    mydisp = @(s) disp(s);
else
    mydisp = @(s) [];
end

S_old=[];
n_it=1;
mydisp('Iteration 1')
[S,Q]=genlouvain(B,limit,verbose,randord,randmove,S0);

mydisp('')

Q_old=-inf;
while ~isequal(S,S_old)
    n_it=n_it+1;
    S_old=S;
    Q_old=Q;
    
    mydisp(sprintf('Iteration %u',n_it))
    if ~isempty(postprocessor)
        S=postprocessor(S);
    end
    [S,Q]=genlouvain(B,limit,verbose,randord,randmove,S);
    mydisp(sprintf('Improvement in modularity: %f\n',Q-Q_old))
end

end
