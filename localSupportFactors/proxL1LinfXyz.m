function proxL1Linf = proxL1LinfXyz(z, lambda)
%
% This method computes the proximal mapping of the l1/l_inf norm.
%
% It is a vectorised implementation of John Duchi's code (see
% https://web.stanford.edu/~jduchi/projects/DuchiShSiCh08.html ). For fixed
% D=3, the runtime complexity is linear in N.
%
% Author: Florian Bernard (2016)
%
    Nall = numel(z);
    N = Nall/3;
    XG = reshape(z,N,3)';
    XGsorted = sort(abs(XG), 'descend');
    XGcumsum = cumsum(XGsorted);

    RHO = ones(1,N);
    logcond3 = XGsorted(3,:) > (XGcumsum(3,:) - lambda)./3;
    RHO(1,logcond3) = 3;
    logcond2 = XGsorted(2,:) > (XGcumsum(2,:) - lambda)./2;
    RHO(1,~logcond3 & logcond2) = 2;

    linIdx = sub2ind(size(XG), RHO, 1:N);
    
    tmp = (XGcumsum(linIdx) - lambda)./RHO;
    THETA = tmp;
    THETA(tmp<0) = 0;

    abszminusTheta = abs(z) - repmat(THETA, 1, 3)';
    posIdx = abszminusTheta>0;
    proxL1Linf = z;
	
	% prox = Id - projOntoL1Ball
    proxL1Linf(posIdx) = z(posIdx) - sign(z(posIdx)).*abszminusTheta(posIdx);     
end