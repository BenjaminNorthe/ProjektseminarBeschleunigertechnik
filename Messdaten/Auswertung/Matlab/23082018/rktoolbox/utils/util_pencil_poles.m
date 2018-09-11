function [xi, c] = util_pencil_poles(K, H, s0)
%UTIL_PENCIL_POLES    Computes the poles of a (block) pencil (H, K).
%
% xi = util_pencil_poles(K, H) returns the poles of the
% (block) pencil (H, K). The optional parameter s0 is the size of the
% first block (s0=1 in the non-block case). If s0>1 and deflation occurred 
% in the pencil, the function util_pencil_blocksizes is used to first infer
% the sizes of all blocks.

if nargin == 2
    s0 = 1;
end

m = size(H, 2); 
xi = H(1:0, 1:1); % initialize to work with mp/vpa

if s0 ~=1
    c = util_pencil_blocksizes(K,H,s0);  
    c = cumsum([1,c]);
    if((c(2)-c(1)) ~= (c(end)-c(end-1)))
        error('Can''t rerun when deflation occurs')
    end    
else
    c = 1:m+2;
end

j = 1;
while c(j) <= m
    if c(j+2) <= m+1
        if double(H(c(j+2):c(j+3)-1, c(j):c(j+1)-1)) == zeros(c(j+3)-c(j+2),c(j+1)-c(j));
            X = H(c(j+1):c(j+2)-1, c(j):c(j+1)-1)./K(c(j+1):c(j+2)-1, c(j):c(j+1)-1);
            X = X(:); 
            l = length(X);
            X = X(setdiff(1:l,find(isnan(X)==1)));
            xi(1,j) = mean(X); % Remove NaN and find mean.
        else
            [HH, KK] = qz(H(c(j+1):c(j+3)-1, c(j):c(j+2)-1), K(c(j+1):c(j+3)-1, c(j):c(j+2)-1));
            cxi = diag(HH)./diag(KK);
            if not(or(isinf(cxi(1)), isreal(cxi(1))));
                cxi_r = mean(real(cxi));
                cxi_i = mean(abs(imag(cxi)));
                % Specifically choose pole with negative imaginary part
                % first. This corrresponds to select function when using
                % ordqz.
                xi(1, j)   = cxi_r - (cxi_i)*1i; 
                xi(1, j+1) = cxi_r + (cxi_i)*1i;
            else
                xi(1,j) = min(cxi);
                xi(1,j+1) = max(cxi);
            end
            j = j+1;
        end
    else
            X = H(c(j+1):c(j+2)-1, c(j):c(j+1)-1)./K(c(j+1):c(j+2)-1, c(j):c(j+1)-1);
            X = X(:); 
            l = length(X);
            X = X(setdiff(1:l,find(isnan(X)==1)));
            % Scaling causes issues with plus and minus inf, mean gave NaN. 
            if isinf(X)
                xi(1,j) = inf;
            else
                xi(1,j) = mean(X); % Remove NaN and find mean.
            end
    end
    j = j+1;
end
end
