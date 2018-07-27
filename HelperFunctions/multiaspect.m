function [B, twom] = multiaspect(A, gamma, omega, type)
    aspects = size(A);
    if numel(omega) == 1
        omega = repmat(omega, numel(aspects), 1);
    end
    if numel(gamma) == 1
        gamma = repmat(gamma, aspects);
    end
    na = prod(aspects);
    [B, twom] = arrayfun(@(A, g) modularity(A{1}, g), A, gamma, ...
                         'UniformOutput', false);
    B = blkdiag(B{:});
    C = zeros(prod(aspects));
    for a = 1:numel(aspects)
        C = C + kron(kron(eye(prod(aspects(a+1:end))),coupling(aspects(a), type(a))*omega(a)),...
                     eye(prod(aspects(1:a-1))));
    end
    C = kron(C, eye(length(A{1})));
    B = B + C;
    twom = sum([twom{:}]) + sum(sum(C));
end

function C = coupling(size, type)
    switch type
        case {'m', 'c'}
            C = ones(size)-eye(size);
        case {'t', 'o'}
            C = diag(ones(size-1,1),1);
            C = C + C';
    end
end
