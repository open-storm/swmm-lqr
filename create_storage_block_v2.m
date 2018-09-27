function [A] = create_storage_block_v2(T,As,lambda)

if nargin < 3 || isnan(lambda)
    lambda = 0;
end


A = (1-lambda*T/As);