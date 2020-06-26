function Dist = computeAxonsSuperpositionMatrix(axon_collection, dims)

amin = zeros(length(axon_collection),2);
amax = zeros(length(axon_collection),2);
for k = 1:length(axon_collection)
    amin(k,:) = min(axon_collection(k).data);
    amax(k,:) = max(axon_collection(k).data);
end

matrix_dims = length(axon_collection);
Dist = zeros(matrix_dims, matrix_dims);
for k = 1:matrix_dims
    for l = k+1:matrix_dims
%         Dist(k,l) = areAxonsSuperposed(axon_collection(k).data, axon_collection(l).data, dims);
        if any(amin(k,:) > amax(l,:)) || any(amax(k,:) < amin(l,:))
            Dist(k,l) = 0;
        else
            Dist(k,l) = areAxonsSuperposed_fast(axon_collection(k).data, axon_collection(l).data, dims); % for speed
        end
%         Dist(l,k) = Dist(k,l);
    end
end
Dist = Dist + Dist.';
end