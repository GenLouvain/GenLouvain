function S = tidy_config(S)
%This function remains almost identical to that originally written by
%Stephen Reid for his greedy.m code.
%   tidy up S i.e.  S = [2 4 2 6] -> S = [1 2 1 3]
T = zeros(size(S));
m = 0;
for i = 1:numel(S)
    if T(i) == 0
        T(S==S(i)) = m + 1;
        m = m + 1;
    end
end
S = T;
