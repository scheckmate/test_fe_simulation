%Based on:
%"Programing the Finite Element Method with Matlab"
%"Jack Chessa"
function [ ] = CalcElastoDynamicsImplicit( numx,numy,delta_t,t_end )
%Calculation of a 2D linear elastodynamic problem
%Calculation of a 2D linear elastodynamic problem at the moment only for a beam
%numx is number of elements in x-direction
%numy is number of elements in y-direction
%delta_t is the timestep
%at t_end the calculation stops
tic;
disp([num2str(toc),' START'])
L = 50; %length of beam (x-direction)
c = 10; %half of height of beam (y-direction)
t = 0.0;
iter = 0;
figure;
ax = axes;
[C,E0,nu0,rho] = GetMaterial();
disp([num2str(toc),' Create Mesh'])
%create mesh of triangles
[Nodes,Elements,~,~] = GetMesh( numx,numy,c,L );
%plot mesh
PlotMesh(Nodes,Elements,c,L,0,ax);

nn = size(Nodes,1);

disp([num2str(toc),' Create Mass Matrix'])
%create stiffness matrix K
MassMatrix = CalcMassMatrix(nn,Nodes,Elements,rho);

disp([num2str(toc),' Create Stiffness Matrix'])
%create stiffness matrix K
StiffnessMatrix = CalcStiffnessMatrix(nn,Nodes,Elements,C);

disp([num2str(toc),' Create Boundary Vector'])
%create vector f from boundary integrals
%Boundary stresses are zero in our example
%BoundaryVector = CalcBoundary(nn,Nodes,BoundaryLeft,BoundaryRight,P,I0,c,L);

[ NodesFixed, NodesMoved, Fixed, Moved ] = GetDirichlet(Nodes,c,L);
uDofsFixed = NodesFixed;
vDofsFixed = NodesFixed + nn;
vDofsMoved = NodesMoved + nn;
%uDofsMoved = NodesMoved;

U_new = zeros(2*nn,1);
U_old = zeros(2*nn,1);
V = zeros(2*nn,1);
A = zeros(2*nn,1);

while(t < t_end)
%while(iter < maxiter)
    t = t + delta_t;
    %Boundary stresses are zero in our example
    RHS = MassMatrix * (4/delta_t^2 * U_new + 4 / delta_t * V + A);
    Matrix = 4/delta_t^2 * MassMatrix + StiffnessMatrix;
    %apply dirichlet boundary conditions (update Matrix and RHS accordingly)
    bcwt = mean(diag(Matrix)); %value needed to "maintain" condition number
    %udofs = NodesFixedX;
    %vdofs = NodesFixedY + nn;
    RHS = RHS - Matrix(:,uDofsFixed) * Fixed;
    RHS = RHS - Matrix(:,vDofsFixed) * Fixed;
    RHS = RHS - Matrix(:,vDofsMoved) * (t * Moved);
    RHS(uDofsFixed) = Fixed * bcwt;
    RHS(vDofsFixed) = Fixed * bcwt;
    RHS(vDofsMoved) = Moved * t * bcwt;
    Matrix(uDofsFixed,:) = 0;
    Matrix(vDofsFixed,:) = 0;
    Matrix(vDofsMoved,:) = 0;
    Matrix(:,vDofsMoved) = 0;
    Matrix(:,uDofsFixed) = 0;
    Matrix(:,vDofsFixed) = 0;
    Matrix(uDofsFixed,uDofsFixed) = bcwt*speye(length(uDofsFixed));
    Matrix(vDofsFixed,vDofsFixed) = bcwt*speye(length(vDofsFixed));
    Matrix(vDofsMoved,vDofsMoved) = bcwt*speye(length(vDofsMoved));
    U_new = Matrix \ RHS;

    A = 4/delta_t^2 * (U_new - U_old) - 4/delta_t * V - A;
    V = 2/delta_t * (U_new - U_old) - V;
    U_old = U_new;
    iter = iter + 1;
    %StiffnessMatrix = CalcStiffnessMatrix(nn,Nodes+[U(1:nn) U(nn+1:2*nn)],Elements,C);
    %StiffnessMatrix(33,50);
end
%Calculate the stresses
disp([num2str(toc),' Calculate Stresses'])
Stresses = CalcStresses(Nodes,Elements,U_new,nn,C);
VonMisesStresses = zeros(size(Elements,1),1);
VonMisesStresses(:,1) = sqrt(Stresses(:,1).^2 + Stresses(:,2).^2 - Stresses(:,1).*Stresses(:,2) + 3*Stresses(:,3).^2);
%plot updated mesh
PlotMesh(Nodes+[U_new(1:nn) U_new(nn+1:2*nn)],Elements,c,L,VonMisesStresses,ax);
disp(['time: ',num2str(t)])
disp(['number of iterations: ', num2str(iter)])
disp(['norm of solution u: ', num2str(norm(U_new))])
disp(['maximum von mises stress: ', num2str(max(VonMisesStresses))])
disp([num2str(toc),' Calculation finished'])
end
