function [ StiffnessMatrix ] = CalcStiffnessMatrix( nn,nodes,elements,C )
%Calculate stiffness matrix

StiffnessMatrix = sparse(2*nn,2*nn);
ninel = size(elements,2); %number of nodes in element (e.g. 3 for triangles)
%get integration points and weights for numeric integration
A = NumIntegration(2,2);
%extract integration points and weights
intpts = A(:,1:2);
intweights = A(:,3);

for e=1:size(elements,1) %loop over elements
   sctr = elements(e,:); %node numbers of element needed for assembling
   sctrB = [sctr sctr+nn];
   
   for q=1:size(intweights) %loop over number of integration points
      pt = intpts(q,:);
      weight = intweights(q);
      [~,dNdxi] = LagrangeBasis(2,pt); %get "local" gradient of shape functions
      J0 = nodes(sctr,:)' * dNdxi; %calculate Jacobi matrix
      detJ0 = det(J0);
      invJ0 = inv(J0);
      dNdx = dNdxi * invJ0; %calculate "global" gradient of shape functions
      
      %calculate B matrix
      B = zeros(3,2*ninel);
      B(1,1:ninel) = dNdx(:,1)';
      B(2,ninel+1:2*ninel) = dNdx(:,2)';
      B(3,1:ninel) = dNdx(:,2)';
      B(3,ninel+1:2*ninel) = dNdx(:,1)';
      
      %calculate and assemble values of stiffness matrix
      StiffnessMatrix(sctrB,sctrB) = StiffnessMatrix(sctrB,sctrB) + B' * C * B * detJ0 * weight;
   end
end
end

