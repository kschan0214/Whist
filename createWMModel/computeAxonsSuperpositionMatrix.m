function Dist = computeAxonsSuperpositionMatrix(axon_collection, dims)

% precompute the amin and amax, KC
amin = zeros(length(axon_collection),2);
amax = zeros(length(axon_collection),2);
for k = 1:length(axon_collection)
    amin(k,:) = min(axon_collection(k).data);
    amax(k,:) = max(axon_collection(k).data);
end
% precompute the logical OR statements, KC
is_amin_gt_amax = zeros(length(axon_collection),length(axon_collection));
for k = 1:length(axon_collection)
    is_amin_gt_amax(:,k) = any(amin(k,:)>amax,2);
end
is_amin_gt_amax = max(is_amin_gt_amax,is_amin_gt_amax.'); % equivalent to is_amin_gt_amax(l,k) || is_amin_gt_amax(k,l)

matrix_dims = length(axon_collection);
Dist = zeros(matrix_dims, matrix_dims);
for k = 1:matrix_dims
    for l = k+1:matrix_dims
%         Dist(k,l) = areAxonsSuperposed(axon_collection(k).data, axon_collection(l).data, dims);
%         if any(amin(k,:) > amax(l,:)) || any(amax(k,:) < amin(l,:))
        if ~is_amin_gt_amax(l,k)
%             Dist(k,l) = 0;
%         else
            Dist(k,l) = areAxonsSuperposed_fast(axon_collection(k).data, axon_collection(l).data, dims); % for speed
        end
%         Dist(l,k) = Dist(k,l);
    end
end
Dist = Dist + Dist.';
end