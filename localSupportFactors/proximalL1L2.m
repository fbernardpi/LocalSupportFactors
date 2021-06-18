function [X,normX,dummyOut] = proximalL1L2(...
    Y,lam,k,~,~,~, ~, ~)
%
% This method computes the proximal operator of the (weighted) sum of the
% l1 norm and l2 norm. Note that in our application the weight for the l1
% norm is set to 0 in order to obtain the norm for A, as described in eq.
% (5) in the paper.
%
% It serves as interface for the FactorTVL1L2. It is based on
% proximalTVL1L2.m, where details on the arguments can be found.
%
% Authors: Ben Haeffele, modified by Florian Bernard (2016)
%

[~,M] = size(Y);

if length(lam)==1
    lam = lam*ones(M,1);
else
    lam = lam(:);
end

if size(k,1)==1
    k = repmat(k,M,1);
end

X = Y; 

if ( any(lam) )
    X = zeros(size(Y));
    for c=1:M
        y = Y(:,c);
        l1Lambda = k(c,1)*lam(c);
        
        X(:,c) = proxL1(y, l1Lambda);
    end
end
temp_norm = [];
for c=1:M
    l1Norm = norm(X(:,c),1);

    temp_norm = [temp_norm, l1Norm];
end

normX = diag(k(:,1)*temp_norm)';

% Now project onto the L2 ball (see Theorem 3)
L2 = sqrt(sum(X.^2,1)); %calculate L2 norms

idx_nz = (L2>0) & (k(:,2)>0)'; %find indexes we need to update

if sum(idx_nz)
    normX(idx_nz) = normX(idx_nz) + k(idx_nz,2)'.*L2(idx_nz); %add L2 norm
	
	% block soft-thresholding (see Parikh's "Proximal Algorithms", Ch. 6.5.1)
    scl = 1-min(L2(idx_nz),lam(idx_nz)'.*k(idx_nz,2)')./L2(idx_nz);

    %scale columns of X to finish projection onto L2 ball
    X(:,idx_nz) = bsxfun(@times,X(:,idx_nz),scl);
    
    normX(idx_nz) = normX(idx_nz).*scl; %update norms
end

dummyOut = nan(0,M);