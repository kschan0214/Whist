
function out = areAxonsSuperposed_alternative(axon1, axon2, dims)

% a1min = axon1.amin;
% a1max = axon1.amax;
% 
% a2min = axon2.amin;
% a2max = axon2.amax;

% a1min = min(axon1);
% a1max = max(axon1);
% 
% a2min = min(axon2);
% a2max = max(axon2);

% if ((sum(a1min > a2max) >= 1) || (sum(a1max < a2min) >= 1))
% if ((sum(axon1.amin > axon2.amax) >= 1) || (sum(axon1.amax < axon2.amin) >= 1))
%     out = 0;
%     return;
% end

ind1 = round(axon1);
sub1 = sub2ind(dims,ind1(:,1),ind1(:,2));  

ind2 = round(axon2);
sub2 = sub2ind(dims,ind2(:,1),ind2(:,2));  

C = intersect(sub1,sub2);
% if (length(C) > 0)
if ~isempty(C)
    out = 1;
else 
    out = 0;
end
end



