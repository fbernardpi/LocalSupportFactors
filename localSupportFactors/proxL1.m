function l1prox = proxL1(z,lambda)
%
% This method computes the proximal mapping of the l1 norm.
%
% Author: Florian Bernard (2016)
%
    minusLambdaIdx = z >= lambda;
    plusLambdaIdx = z <= -lambda;
    l1prox = zeros(size(z));
    l1prox(minusLambdaIdx) = z(minusLambdaIdx) - lambda;
    l1prox(plusLambdaIdx) = z(plusLambdaIdx) + lambda;
end