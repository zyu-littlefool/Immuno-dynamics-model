function [H, unique_gene] = get_Shannon_2(X, xE, xE_2)

    unique_gene = struct('seednum', 1, 'count', 0, 'frequen', 0, ...
                         'E_chenge_wt', 0, 'E_chenge_var', 0, 'char', "0");

    if isempty(X)
        H = 0; return;
    end

    if ~islogical(X)
        X = logical(X); 
    end

    [uniqMat, ~, idx] = unique(X, 'rows', 'stable');
    
    nPattern = size(uniqMat, 1);
    counts = accumarray(idx(:), 1, [nPattern, 1]);
    freqs = counts / numel(idx);

    H = -sum(freqs .* log(freqs + eps));

    unique_gene.E_chenge_wt = double(uniqMat) * xE;
    unique_gene.E_chenge_var = double(uniqMat) * xE_2;

    unique_gene.char = cellstr(char(double(uniqMat) + 48));

    unique_gene.seednum = nPattern;
    unique_gene.count = counts;
    unique_gene.frequen = freqs;
end