function [ Stresses ] = CalcStresses( Nodes, Elements, U, nn, C )
%Calculate stresses from result vector U

Stresses = zeros(size(Elements,1),3);

for e=1:size(Elements,1) %loop over elements
    sctr = Elements(e,:); %node numbers of element
    sctrB = [sctr sctr+nn];
    ninel = length(sctr); %number of nodes in element
    stresspoint = [1/3 1/3]; %point in reference triangle at which the stresses are calculated
    
    [~,dNdxi] = LagrangeBasis(2,stresspoint); %get "local" gradient of shape functions
    J0 = Nodes(sctr,:)' * dNdxi; %calculate Jacobi matrix
    invJ0 = inv(J0);
    dNdx = dNdxi * invJ0; %calculate "global" gradient of shape functions
    
    %calculate B matrix
    B = zeros(3,2*ninel);
    B(1,1:ninel) = dNdx(:,1)';
    B(2,ninel+1:2*ninel) = dNdx(:,2)';
    B(3,1:ninel) = dNdx(:,2)';
    B(3,ninel+1:2*ninel) = dNdx(:,1)';
    
    strain=B*U(sctrB);
    Stresses(e,:) = C * strain;
end

end

