function [Phi, A, objective] = computeLocalSupportDeformationModel(X, params)
%
% This method interfaces the structured low-rank matrix factorisation
% framework contained in the folder FactorTVL1L1.
%
% Input: X                   (3*N-by-K matrix) shape vertices that are in
%                            correspondence
%        params              Matlab struct that specificies various
%                            parameters. Parameters marked with (*) are
%                            obligatory.
%    (*) params.rank         (scalar) Maximum rank of the solution
%    (*) params.C            (N-by-N matrix) vertex connectivity, where the
%                            element C(i,j) denotes the affinity of vertex i
%                            and j. A high value indicates that the deformation
%                            of vertex i and j should be similar, whereas 0
%                            indicates that there is no dependency.
%        params.lambda       (scalar) Weight for the entire regularisation 
%                            term
%                            (default 64)
%        params.lambdaG      (scalar) Weight for the graph-based l2 norm 
%                            term of Phi
%                            (default 1)
%        params.lambdainf    (scalar) Weight for the l1/l_inf norm term of
%                            Phi
%                            (default 2)
%        params.lambda1      (scalar) Weight for the l1 norm term of Phi
%                            (default 1)
%        params.lambda2      (scalar) Wweight for the l2 norm term of Phi
%                            (default 1)
%        params.lambdaA      (scalar) Weight for the l2 norm term of A
%                            (default 1e-4)
%        params.display      (0 or 1) toggle visualisation
%                            (default 1)
%        params.verbose      (0 or 1) toggle textual output
%                            (default 1)
%        params.nIterations  (scalar) number of iterations
%                            (default 100)
%
% Author: Florian Bernard (2016)
%

% parse params
M = params.rank;
C = params.C;

if ( ~isfield(params, 'lambda') )
	params.lambda = 64;
end
if ( ~isfield(params, 'lambdaG') )
	params.lambdaG = 1;
end
if ( ~isfield(params, 'lambdainf') )
	params.lambdainf = 2;
end
if ( ~isfield(params, 'lambda1') )
	params.lambda1 = 1;
end
if ( ~isfield(params, 'lambda2') )
	params.lambda2 = 1;
end
if ( ~isfield(params, 'lambdaA') )
	params.lambdaA = 1e-4;
end
if ( ~isfield(params, 'nIterations') )
	params.nIterations = 100;
end
if ( ~isfield(params, 'display') )
	params.display = 1;
end

% set proximal operators
params.proximalOperatorA = @proximalStructuredSparsityDualForwardBackward;
params.proximalOperatorZ = @proximalL1L2;

% enable warmstart, which improves speed of proximal operator computation
params.save_dual = true;

K = size(X,2);
Nall = size(X,1);
N = Nall/3;

% normalise X
meanShape = mean(X,2);

Xcentred = bsxfun(@minus, X, meanShape);
scaling = std(Xcentred(:));
X = Xcentred./scaling;

% remove symmetry in graph to reduce runtime. Also, ignore diagonal
if ( norm(C-C', 'fro') < 1e-8 )
	C = 2*triu(C,1);
end

% setup weights from graph matrix C
[v1,v2] = find(C);
graphNeighbourGroup = [v1,v2];
weights = full(C(sub2ind(size(C),v1,v2)));

graphNeighbourGroup = [graphNeighbourGroup; ...
	graphNeighbourGroup+N; ...
	graphNeighbourGroup+2*N];
weights = repmat(weights,3,1);

params.idx = {graphNeighbourGroup};
params.idx_indexA = 1;
params.idx_indexZ = [];
params.idx_indexB = [];

params.weights = {sqrt(weights)};

% normalise parameters to make them problem-size invariant
params.lam = params.lambda*Nall*K/M;
params.kA = [params.lambdaG/sqrt(numel(weights)) ...
	params.lambdainf/sqrt(N) ...
	params.lambda1/sqrt(Nall) ...
	params.lambda2/sqrt(Nall)];
params.kZ = [0 params.lambdaA/sqrt(K)];


if ( params.display )
	figure;
	colorbar;
	caxis auto;
	
	titleStr = ['$\lambda = ' num2str(params.lambda) ...
		'; \lambda_A = ' num2str(params.lambdaA) ...
		'; \lambda_1 = ' num2str(params.lambda1) ...
		'; \lambda_2 = ' num2str(params.lambda2) ...
		'; \lambda_{\infty} = ' num2str(params.lambdainf) ...
		'; \lambda_{\mathcal{G}} = ' num2str(params.lambdaG) '$'];
	title(titleStr, 'Interpreter','latex', 'FontSize', 14);
	
	hold on;
end

% run the structured low-rank matrix factorisation framework
[Phi,A,~,objective,~,~,~,L_stats] = ...
	FactorTVL1L2_v1(X, [],[],params, params.nIterations);

% use l2 norm for selecting which columns are small and thus
% set to 0
A_norm = sqrt(sum(Phi.*Phi,1));
Z_norm = sqrt(sum(A.*A,1));

% first, remove factors that are close to 0
si = 1:size(A,2);
norms = A_norm.*Z_norm;
si = si(norms >= 1e-8);

Phi = Phi(:,si);
A = A(:,si);

% normalise such that weights have std 1
ZStd = std(A',0,2);
nonzeroidx = logical(ZStd);
A(:,nonzeroidx) = bsxfun(@times, A(:,nonzeroidx)', 1./ZStd(nonzeroidx))';
Phi(:,nonzeroidx) = bsxfun(@times, Phi(:,nonzeroidx), ZStd(nonzeroidx)');

Phi = Phi*scaling;
end