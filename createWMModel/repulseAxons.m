function [axon_collection, FVF_current] = repulseAxons(axon_collection, FVF_expected, tol, mask, plot_model)

[~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);

if FVF_current < FVF_expected
    error('The current FVF before axon dispersions is lower than the expected FVF and cannot be reach')
end

% Define model boundaries 
maxSize = 0;

N = length(axon_collection);
for k = 1:N
    currentSize = max(max(axon_collection(k).data)  - min(axon_collection(k).data ));
    if (currentSize > maxSize)
        maxSize = currentSize;
    end
end

dims = size(mask);
[mask_x_index, mask_y_index] = find(mask);
mask_min_x = min(mask_x_index);
mask_max_x = max(mask_x_index);
mask_min_y = min(mask_y_index);
mask_max_y = max(mask_y_index);

edge = 2;
remove_edge_min_x = max(mask_min_x - (maxSize + edge), (maxSize + edge));
remove_edge_max_x = min(mask_max_x + (maxSize + edge), dims(1) - (maxSize + edge));
remove_edge_min_y = max(mask_min_y - (maxSize + edge), (maxSize + edge));
remove_edge_max_y = min(mask_max_y + (maxSize + edge), dims(2) - (maxSize + edge));

% get axon position
pts = cat(1,axon_collection(:).Centroid);
dims = size(mask);

% step size of moving axons
step = 0.01;

iter_total = 0;
max_iter = 1e6; % set maximimum iteration to avoid being trapped

% while or((FVF_current < FVF_expected - tol),(FVF_current > FVF_expected))
while and(abs(FVF_current - FVF_expected) > tol, iter_total < max_iter)
    iter_total = iter_total + 1;

    display(['currentFVF : ' num2str(FVF_current)]);

    axon_collection_old = axon_collection;
    pts_old             = pts; % KC - bug fix: pts also need to be restored array size not the same
    pts_directions = (pts - [dims(1)/2 dims(2)/2]) ;
    pts = pts + step*pts_directions;
    
    k = 1;
    while k < length(axon_collection)
        axon_collection(k).data = axon_collection(k).data - axon_collection(k).Centroid + pts(k,:);
        
        % Remove axons out of boundaries    
        if ((min(pts(k, 1)) < remove_edge_min_x) || (max(pts(k, 1)) > remove_edge_max_x) || ...
                (min(pts(k, 2)) < remove_edge_min_y) || (max(pts(k, 2)) > remove_edge_max_y))
            
            axon_collection(k) = [];
            pts(k,:) = [];
        else
            % moving axons
            axon_collection(k).Centroid = pts(k,:);
            k = k + 1;
        end      
    end
    
    % get current FVF after the last operation
    if (mod(iter_total, 5) == 0)
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, plot_model);
%         pause(0.01);
        drawnow
    else
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, 0);
    end
    
    % lower bound condition
    % restore iteration with half the step size if exceeds the lb
    if (FVF_current < FVF_expected - tol)
        axon_collection = axon_collection_old;
        pts             = pts_old; % KC - bug fix: pts also need to be restored array size not the same
        [~, ~, FVF_current] = createModelFromData(axon_collection, mask, 0);
        step = step/2;
    end  
end

if iter_total == max_iter
    warning('Maximum no. of iteration reached when repulsing axons.');
end

end