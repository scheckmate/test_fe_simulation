function [ MassMatrix ] = CalcMassMatrix( nn,nodes,elements,rho )
%Calculate stiffness matrix

MassMatrix = sparse(2*nn,2*nn);
%get integration points and weights for numeric integration
A = NumIntegration(2,2);
%extract integration points and weights
intpts = A(:,1:2);
intweights = A(:,3);

for e=1:size(elements,1) %loop over elements
   sctr = elements(e,:); %node numbers of element needed for assembling
   
   for q=1:size(intweights) %loop over number of integration points
      pt = intpts(q,:);
      weight = intweights(q);
      [N,dNdxi] = LagrangeBasis(2,pt); %get "local" gradient of shape functions
      J0 = nodes(sctr,:)' * dNdxi; %calculate Jacobi matrix
      detJ0 = det(J0);
      
      mQPt = rho * (N * N') * detJ0 * weight;
      %calculate and assemble values of mass matrix
      MassMatrix(sctr,sctr) = MassMatrix(sctr,sctr) + mQPt;
      MassMatrix(sctr+nn,sctr+nn) = MassMatrix(sctr+nn,sctr+nn) + mQPt;
   end
end
end
