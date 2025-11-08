function [L2, L1, linf] = norms(v)

    L2 = sqrt(sum(v.^2));
    L1 = sum(abs(v));
    linf = max(abs(v));
end