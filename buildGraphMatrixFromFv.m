function G = buildGraphMatrixFromFv(fvOrFaces)
%
% This method creates a binary or weighted graph, depending on the input.
% If fvOrFaces is a struct with the fields "vertices" and "faces", then the
% euclidean distance between neighbouring vertices is used as graph weight.
% If fvOrFaces is only a faces variable, the output is a binary graph of
% the connectivity specified by faces.
%
% Author: Florian Bernard (2016)
%

if ( isstruct(fvOrFaces) )
	vertexFeatures = fvOrFaces.vertices;
	faces = fvOrFaces.faces;
	weightedGraph = true;
else
	faces = fvOrFaces;
	weightedGraph = false;
end

N = numel(unique(faces(:)));

G = sparse(N,N);
for f=1:size(faces,1)
	neighs = faces(f,:);

	% for all edges do
	edges = nchoosek(neighs,2);

	for e=1:size(edges,1)
		v1Idx = edges(e,1);
		v2Idx = edges(e,2);

		if ( weightedGraph )
			v1 = vertexFeatures(v1Idx,:);
			v2 = vertexFeatures(v2Idx,:);

			vertDist = norm(v1-v2);

			G(v1Idx,v2Idx) = vertDist;
			G(v2Idx,v1Idx) = vertDist;
		else
			G(v1Idx,v2Idx) = 1;
			G(v2Idx,v1Idx) = 1;
		end
	end
end
end