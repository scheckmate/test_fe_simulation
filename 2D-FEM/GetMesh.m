function [ Nodes, Elements, BoundaryLeft, BoundaryRight ] = GetMesh( numx, numy, c, L )
%Calculate a 2D mesh of triangles of a beam
%numx is number of elements in x-direction
%numy is number of elements in y-direction
%2*c is height of beam
%L is length of beam
%The coordinates of the nodes are returned in the vector Nodes
%The number of the nodes for each element is returned in the vector
%Elements
%The number of the nodes for the boundary elements is returned in the
%vector BoundaryLeft (x=0) and Boundary Right (x=L) respectively

nodesx = numx +1;
nodesy = numy +1;
Nodes = zeros(nodesx*nodesy,2);
%i=0:numy;
j=0:numx;
x_coords = j*L/numx;
Nodes(:,1) = repmat(x_coords,1,nodesy);
for m=0:numy
Nodes(m*nodesx+1:(m+1)*nodesx,2) = c-m*2*c/numy;
%Nodes(i.*nodesx+1:(i+1).*nodesx,2) = c-i*2*c/numy;
end
Elements = zeros(2*numx*numy,3);
for l=0:numy-1
    for k=1:numx
        Elements(l*2*numx+2*k-1,:) = [k+l*nodesx k+(l+1)*nodesx k+1+l*nodesx];
        Elements(l*2*numx+2*k,:) = [k+(l+1)*nodesx k+1+(l+1)*nodesx k+1+l*nodesx];
    end
end

BoundaryLeft = zeros(numy,2);
BoundaryRight = zeros(numy,2);

for i=1:numy
    BoundaryLeft(i,:) = [(i-1)*nodesx+1 i*nodesx+1];
    BoundaryRight(i,:) = [(numy-i+2)*nodesx (numy-i+1)*nodesx];
end

end

