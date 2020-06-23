function Dist = computeAxonsSuperpositionMatrix(axon_collection, dims)

for k = 1:length(axon_collection)
    axon_collection(k).amin = min(axon_collection(k).data);
    axon_collection(k).amax = max(axon_collection(k).data);
end

matrix_dims = length(axon_collection);
Dist = zeros(matrix_dims, matrix_dims);
for k = 1:matrix_dims
    for l = k+1:matrix_dims
%         Dist(k,l) = areAxonsSuperposed(axon_collection(k).data, axon_collection(l).data, dims);
        Dist(k,l) = areAxonsSuperposed_alternative(axon_collection(k), axon_collection(l), dims); % for speed
        Dist(l,k) = Dist(k,l);
    end
end
end