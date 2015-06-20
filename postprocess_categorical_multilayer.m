function S=postprocess_categorical_multilayer(S,T,max_coms)

% post-process an NxT multilayer partition S with N nodes and T layers to
% maximise multiplex persistence without changing the community structure 
% within layers. The algorithm iterates through layers in random order and 
% uses the Hungarian algorithm to solve the optimal assignment
% problem for each layer. The algorithm stops once it cannot improve the 
% assignment for any layer. Note that this procedure always increases multiplex
% multilayer modularity for any non-zereo value of omega.
%
% The output partition is stochastic due to random order and unlike in the
% temporal case it is not guaranteed to have optimal persistence for the
% given intra-layer community assignments.

if nargin<2||isempty(T)
    T=size(S,2);
end
N=numel(S)/T;


if nargin<3||isempty(max_coms)
    max_coms=inf;
end

[N0,T0]=size(S);

if max(S(:))<max_coms % don't do anything if to many communities for performance

% tidy assignment
[~,~,S]=unique(S);    
S=reshape(S,N,T);
S_new=S;

p0=0;
p1=multiplex_persistence(S); % checking

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
        if S2(j)~=0
            S_new(Gi(:,j),i)=uc(S2(j)); % update assignment
        else
            S_new(Gi(:,j),i)=max_com+1; % get new label if not assigned
            max_com=max_com+1;
        end
        end
    end
    
    S_new=tidy_config(S_new);
    
    p1=multiplex_persistence(S_new);
    
    fprintf('improvement found: %g\n',p1-p0);
end

% return in original format
S=reshape(S,N0,T0);

end
