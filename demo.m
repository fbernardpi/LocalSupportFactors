%% load and visualise toy data
addpath(genpath(pwd)); % add subfolders to path path

load shapes.mat

vis = 1; % do visualisation

K = numel(shapes); % number of shapes
X = nan(numel(shapes{1}), K); % 3N-by-K matrix of vertex coordinates

if ( vis )
	figure;
	rotate3d;
	m = ceil(sqrt(K));
end

for k=1:K
	currVerts = shapes{k};
	X(:,k) = currVerts(:);
	
	if ( vis )
		fv.vertices = currVerts;
		fv.faces = faces;
	
		subplot(m,m,k);
		hold on;
		patch(fv, 'FaceColor', 'g', 'FaceAlpha', 0.5);
	
		axis equal;
		axis tight;	
	end
end

%% obtain sparse localised deformations
Nall = size(X,1);
N = Nall/3;

% compute mean shape in order to compute the connectivity graph
meanShape = mean(X,2);
meanShape3 = reshape(meanShape,N,3);
meanShape3 = meanShape3./std(meanShape3(:)); 

% construct connectivity graph from normalised mean shape
meanFv.vertices = meanShape3;
meanFv.faces = faces;
G = buildGraphMatrixFromFv(meanFv);

% use Gaussian kernel of (normalised) euclidean distance between
% neighbouring vertices as weights for the graph
meanLmDistance = mean(G(logical(G)));
           
[e1, e2] = find(G);
linIdx = sub2ind(size(G), e1, e2);
Dexp = sparse(size(G,1),size(G,2));
Dexp(linIdx) = exp(-(G(linIdx)./meanLmDistance).^2);

% desired (maximum) rank of solution. Note that factors close to 0 are
% removed, so the final number of factors may be smaller than this number
params.rank = 22;
			
% set weighted graph
params.C = Dexp;
% 
% params.display = 1;
% params.verbose = 1;


% perform the actual fitting
randseed(1);
tic
[factors,weights,obj] = computeLocalSupportDeformationModel(X, params);
toc;

% ... at this place we would run the factor splitting in order to split each
% factor that has more than one active region into multiple factors ...

% sort factors according to l2 norm
factorNorms = sqrt(sum(factors.*factors,1));
[~,idx] = sort(factorNorms, 'descend');
factors = factors(:,idx);
weights = weights(:,idx);

M = size(factors,2);


%% show results
figure;
rotate3d;
m = ceil(sqrt(M));

for i=1:M
	subplot(m,m,i);
	factor_i = reshape(factors(:,i), N, 3);
	
	factorMag = sqrt(sum(factor_i.*factor_i,2));

	lmDisplacementMag = (1/(max(factorMag)-min(factorMag)))*(factorMag-min(factorMag));
	nColors = 1024;
	lmColorIdx = ceil((lmDisplacementMag*(nColors-1))+1);
	cmap = colormap(parula(nColors));
	
	fv.facevertexcdata = cmap(lmColorIdx,:);
	fv.vertices = meanShape3;
	fv.faces = faces;
	
	patch(fv, 'FaceColor', 'interp', 'FaceLighting', 'gouraud',...
		'EdgeColor', 'none', 'EdgeLighting', 'gouraud',...
		'MarkerFaceColor', 'auto', ...
		'MarkerEdgeColor', 'auto');
end
