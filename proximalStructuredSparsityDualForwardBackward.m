function [X,normX,Vdual] = proximalStructuredSparsityDualForwardBackward(...
    Y,lam,k,idx,~,Vdual_init,weights, verbose)
%
% This method computes the proximal operator of the (weighted) sum of the
% l1 norm, l2 norm, l1/l_inf norm and the affine composition of the l2
% norm. See eq. (6) in the paper for the precise definition.
%
% It serves as interface for the FactorTVL1L2. It is based on
% proximalTVL1L2.m, where details on the arguments can be found.
%
% Author: Florian Bernard (2016)
%
[Nall,M] = size(Y);

if length(lam)==1
    lam = lam*ones(M,1);
else
    lam = lam(:);
end

if size(k,1)==1
    k = repmat(k,M,1);
end

saveDual = false;
if nargout>=3
    Vdual = zeros(size(idx,1),M);
    saveDual = true;
end
if ~exist('weights', 'var') 
	weights = ones(size(idx,1),1);
end


X = Y;
idx = double(idx);
nEdges = size(idx,1);
L = sparse(1:nEdges, idx(:,1), weights, nEdges,Nall) - ...
    sparse(1:nEdges, idx(:,2), weights, nEdges,Nall);

N = Nall/3;
if ( any(k(:,2)) )
    kronG = kron(ones(1,3),speye(N));
else
    kronG = 0;
end

if ( any(lam) && any(any(k(:,1:3))) ) 
    Lfro = norm(L, 'fro');
    L = L./Lfro;

    ticTmp = tic;
    for c=1:M
        y = Y(:,c);
        xi = y;

        l2AffineLambda = lam(c)*k(c,1);
        l2AffineLambda = l2AffineLambda*Lfro;
        l1linfLambda = lam(c)*k(c,2);
        l1Lambda = lam(c)*k(c,3);
        l2Lambda = lam(c)*k(c,4);

        %% DUAL FORWARD BACKWARD SCHEME
        maxIt = 20;
        epsilon = 1e-4;
        if ( l2AffineLambda > 0 ) % dual forward backward method required

%% DUALIZATION
            eps = 1-1e-4;
            gamma_k = 2-eps;
            lambda_k = (1+eps)/2;

            conv = [];
            if ( ~isempty(Vdual_init) && any(Vdual_init(:)) ) % non-trivial Vdual_init given
                v_k = Vdual_init(:,c);
            else
                v_k = 0.1*(rand(size(L,1),1)-0.5);
            end
            v_kPrev = inf(size(L,1),1);
            x_k = zeros(size(y));
            it = 0;

            Ladj = L';

            while ( it < maxIt )
                proxfarg = y - Ladj*v_k;
                % compute l1/linf + l1 prox
                %  compute l1 prox first ...
                if ( l1Lambda > 0 )
                    l1prox = proxL1(proxfarg, l1Lambda);
                else
                    l1prox = proxfarg;
                end
                
                % ... then compute l1/linf prox
                if ( l1linfLambda > 0 )
                    x_k = proxL1LinfXyz(l1prox, l1linfLambda);
                else
                    x_k = l1prox;
                end

                if ( l2Lambda > 0 )
                    x_k = proxL2(x_k, l2Lambda);
				end
     
				% update v_k
				Ex = L*x_k;
				proxl2tmp = proxL2(v_k./gamma_k+Ex, l2AffineLambda/gamma_k);
				v_k = v_k + lambda_k*gamma_k*(Ex - proxl2tmp);
	
%                 % using idx/weights instead of L is not faster
%                 proxGarg = v_k + gamma_k*L*x_k;
% 
%                 % direct solution of prox of l2 norm
%                 proxL2NormArg = proxGarg./gamma_k;
%                 proxL2NormGamma = l2AffineLambda/gamma_k;
%                 proxl2tmp = proxL2(proxL2NormArg, proxL2NormGamma);
%                 
%                 proxl2 = proxGarg - gamma_k*proxl2tmp;
% 
%                 v_k = v_k + lambda_k * (proxl2 - v_k);

                if (it > 0 && (norm(v_k-v_kPrev,2) <= epsilon )) 
                    break;
                end
                v_kPrev = v_k;

                it = it+1;
            end

            xi = x_k;
            if ( saveDual )
                Vdual(:,c) = v_k;
            end
            hold on, plot(conv);
        else % no iterative scheme necessary
            if ( l1Lambda > 0 )
                xi = proxL1(y, l1Lambda);
            else
                xi = y;
            end
            
			if ( l1linfLambda > 0 )
				xi = proxL1LinfXyz(xi, l1linfLambda);
			end
			
			if ( l2Lambda > 0 )
				xi = proxL2(xi, l2Lambda);
			end
		end
        X(:,c) = xi;     
    end
    s=toc(ticTmp);
    if ( verbose )
        disp(['prox computed in ' num2str(s) 's']);
    end
end
temp_norm = [];
for c=1:M
    l2AffineNorm = norm(L*X(:,c),2);
    l1l2Norm = sum(sqrt(kronG*(X(:,c).*X(:,c))));
    l1Norm = norm(X(:,c),1);
    l2Norm = norm(X(:,c),2);
    temp_norm = [temp_norm, ...
        [l2AffineNorm; l1l2Norm; l1Norm; l2Norm]];
end

if ( isempty(temp_norm) )
    normX = zeros(M,2);
else
    normX = diag(k(:,1:4)*temp_norm)';
end
