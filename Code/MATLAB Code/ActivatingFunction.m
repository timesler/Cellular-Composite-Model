function [AF] = ActivatingFunction(tree,Ve)

p = NTESparams('double');

AF = zeros(size(Ve));

for i = 1:length(tree.X)
    % Get potential at current node
    Ve_n = Ve(:,i,:);
%     Ri_n = tree.Ri(i);
    
    % Get potential at parent and children nodes
    Ve_adj = Ve(:,logical(tree.dA(i,:)) | logical(tree.dA(:,i)'),:);
    Ri_adj = tree.Ri(logical(tree.dA(i,:)) | logical(tree.dA(:,i)'));
    D_adj = tree.D(logical(tree.dA(i,:)) | logical(tree.dA(:,i)'));
    
    % Get axial resistance (different for soma)
    if tree.R(i) ~= 1
        Ri_n = tree.Ri(i);
    end
    
    for j = 1:length(Ri_adj)
        % if we are at the soma, assume a tapered geometry
        if tree.R(i) == 1
            Ri_n = p.rho_i*tree.seglen(i)./(pi*(tree.D(i)+D_adj(j)).^2/8);
        end
        AF(:,i,:) = AF(:,i,:) + (Ve_adj(:,j,:) - Ve_n)/((Ri_adj(j) + Ri_n)/2);
    end
end