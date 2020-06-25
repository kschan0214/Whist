% function axon_index = myelin2axon_tmp(myelin_index)
function axonCollection = myelin2axon_fast(axonCollection)

nx = 64;
myelin_map_all = zeros(nx,nx*length(axonCollection));

for k= 1:length(axonCollection)     
    
    ind_start   = 1 + (k-1)*nx;

    myelin_index = round(axonCollection(k).data);
    
    indminx = min(myelin_index(:,1))-1;
    indminy = min(myelin_index(:,2))-1;
    
    sizeax = [max(myelin_index(:,1)) - indminx, max(myelin_index(:,2)) - indminy];
    % Initialize myelin with small size
    myelin_map = zeros(sizeax);
    ind = sub2ind(sizeax, myelin_index(:,1) - indminx, myelin_index(:,2) - indminy);
    myelin_map(ind) = 1;
    
    myelin_map_all(1:sizeax(1),ind_start:ind_start+sizeax(2)-1) = myelin_map;
    
    axonCollection(k).myelin = myelin_index;
end
    
axon_map_all = imfill(myelin_map_all, 'holes') - myelin_map_all;

for k= 1:length(axonCollection) 
    ind_start   = 1 + (k-1)*nx;
    
    myelin_index = axonCollection(k).myelin;
    
    indminx = min(myelin_index(:,1))-1;
    indminy = min(myelin_index(:,2))-1;
    
    sizeax = [max(myelin_index(:,1)) - indminx, max(myelin_index(:,2)) - indminy];

    curr_axon = axon_map_all(1:sizeax(1),ind_start:ind_start+sizeax(2)-1);
    
    [xa,ya] = find(curr_axon);
    
    if isempty([xa,ya])
        xa=1;
        ya=1;
    end
    
    axon_index = zeros(length(xa),2); % bug fix: KC 20200625
    axon_index(:,1) = xa + indminx;
    axon_index(:,2) = ya + indminy;
    axonCollection(k).axon = axon_index;
end


% % Get axon region (BG region inside myelin mask)
% axon_map = imfill(myelin_map, 'holes') - myelin_map;
% [xa,ya] = find(axon_map);
% 
% if isempty([xa,ya])
%     xa=1;
%     ya=1;
% end
% 
% % axon_index = zeros(length(xa + indminx),2);
% axon_index = zeros(length(xa),2); % bug fix: KC 20200625
% axon_index(:,1) = xa + indminx;
% axon_index(:,2) = ya + indminy;

end