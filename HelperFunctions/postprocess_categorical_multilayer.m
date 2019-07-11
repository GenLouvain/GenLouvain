function S=postprocess_categorical_multilayer(S,T,max_coms,verbose)
% POSTPROCESS_CATEGORICAL_MULTILAYER post-process an unordered multilayer partition
% to improve persistence when using uniform categorical coupling
% with the multilayer quality function in Mucha et al. 2010.
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
% Call as:
%
%     S=postprocess_categorical_multilayer(S)
%
%     S=postprocess_categorical_multilayer(S,T)
%
%     S=postprocess_categorical_multilayer(S,T,max_coms)
%
% Input:
%
%     S: multilayer partition
%
%     T: number of layers of S (defaults to 'size(S,2)')
%
%     max_coms: only run function when input partition has less than
%         'max_coms' communities, otherwise return input partition.
%         (defauts to 'inf')
%
% Output:
%
%     S: post-processed multilayer partition
%
% The algorithm iterates through layers in a random order and
% uses the Hungarian algorithm to solve the optimal assignment
% problem for each layer. The algorithm stops once it cannot improve the
% assignment for any layer. Note that this procedure always increases
% multilayer modularity with uniform categorical interlayer coupling for
% any non-zereo value of coupling strength omega.
%
% The output partition is stochastic due to random order and unlike in the
% ordinal case it is not guaranteed to have optimal persistence for the
% given intralayer community assignments.
%
%   References:
%
%     Mucha, Peter J., Thomas Richardson, Kevin Macon, Mason A. Porter, and
%     Jukka-Pekka Onnela. "Community Structure in Time-Dependent,
%     Multiscale, and Multiplex Networks," Science 328, 876-878 (2010).
%
%     Bazzi, Marya, Mason A. Porter, Stacy Williams, Mark McDonald, Daniel
%     J. Fenn, and Sam D. Howison. "Community Detection in Temporal
%     Multilayer Networks, with an Application to Correlation Networks",
%     MMS: A SIAM Interdisciplinary Journal 14, 1-41 (2016).


if nargin<2||isempty(T)
    T=size(S,2);
end
N=numel(S)/T;


if nargin<3||isempty(max_coms)
    max_coms=inf;
end

if nargin<4||isempty(verbose)
    verbose=false;
end

if verbose
    mydisp=@(s) disp(s);
else
    mydisp=@(s) [];
end

[N0,T0]=size(S);

if max(S(:))<max_coms % don't do anything if too many communities for performance

% tidy assignment
[~,~,S]=unique(S);
S=reshape(S,N,T);
S_new=S;

p0=-inf;
p1=categorical_persistence(S); % checking

while (p1-p0)>0
    p0=p1;
    % only update if improvement found (don't update in degenerate case
    % where change in persistence is 0)
    S=S_new;
    max_com=max(S_new(:));
    order=randperm(T);
    for i=order
        [ui,~,ei]=unique(S_new(:,i));
        Gi=sparse(1:length(ei),ei,true);
        c=1:T;
        c(i)=[];
        [uc,~,ec]=unique(S_new(:,c));
        ec=reshape(ec,N,T-1);
        overlap=zeros(length(ui),length(uc));
        for j=1:length(ui)
            ecj=ec(Gi(:,j),:);
            Gc=sparse(1:numel(ecj),ecj(:),1,numel(ecj),length(uc));
            overlap(j,:)=full(sum(Gc,1));
        end
        dist=sum(overlap(:))-overlap;
        S2=assignmentoptimal(dist);

        for j=1:length(ui)
        if S2(j)~=0&&overlap(j,S2(j))>0
            S_new(Gi(:,j),i)=uc(S2(j)); % update assignment
        else
            S_new(Gi(:,j),i)=max_com+1; % get new label if not assigned
            max_com=max_com+1;
        end
        end
    end
    p1=categorical_persistence(S_new);

    mydisp(sprintf('improvement found: %g',p1-p0));
end

% return in original format
S=reshape(S,N0,T0);

end
