function [ NodesFixed, NodesMoved, Fixed, Moved ] = GetDirichlet(Nodes,c,L)
%return Dirichlet nodes and the corresponding values
vel = -1.0;

NodesFixed = find(Nodes(:,1) == 0.0);
[~,NodesMoved]=ismember([L,c],Nodes,'rows');
Fixed = zeros(size(NodesFixed,1),1);
%NodesMoved = find(Nodes(:,1) == L);
Moved = vel * ones(size(NodesMoved,1),1);

end
