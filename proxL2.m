function proxL2 = proxL2(z, lambda)
%
% This method computes the proximal mapping of the l2 norm.
%
% Author: Florian Bernard (2016)
%
    normVal = sqrt(sum(z.^2));
    if ( normVal >= lambda )
        proxL2 = (1-lambda/normVal).*z;
    else
        proxL2 = zeros(size(z));
    end
end