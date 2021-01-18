%Based on:
%"Programing the Finite Element Method with Matlab"
%"Jack Chessa"
function [ ] = CalcElastoDynamicsExplicit( numx,numy,alpha,t_end )
%Calculation of a 2D linear elastodynamic problem
%Calculation of a 2D linear elastodynamic problem at the moment only for a beam
%numx is number of elements in x-direction
%numy is number of elements in y-direction
%alpha is a factor with which the critical time step is multiplied (could
%be set to e.g. 0.9)
%at t_end the calculation stops
tic;
disp([num2str(toc),' START'])
L = 50; %length of beam (x-direction)
c = 10; %half of height of beam (y-direction)
t = 0.0;
iter = 0;
figure;
ax = axes;
lumping = 1;
[C,E0,nu0,rho] = GetMaterial();
disp([num2str(toc),' Create Mesh'])
%create mesh of triangles
[Nodes,Elements,~,~] = GetMesh( numx,numy,c,L );
trig_a = L/numx;
trig_b = 2*c/numy;
char_length = 2 * (trig_a * trig_b / 2) / sqrt(trig_a^2 + trig_b^2);
wavespeed = sqrt(E0/(rho*(1-nu0^2)));
delta_t_crit = char_length / wavespeed
delta_t = alpha * delta_t_crit;
%plot mesh
PlotMesh(Nodes,Elements,c,L,0,ax);

nn = size(Nodes,1);

disp([num2str(toc),' Create Mass Matrix'])
%create stiffness matrix K
MassMatrix = CalcMassMatrix(nn,Nodes,Elements,rho);
if(lumping == 1)
    LumpedMass = sum(MassMatrix,2);
end
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

U = zeros(2*nn,1);
V = zeros(2*nn,1);
A = zeros(2*nn,1);

while(t < t_end)
%while(iter < maxiter)
    %Boundary stresses are zero in our example
    %A = MassMatrix \ (BoundaryVector - StiffnessMatrix * U);
    RHS = -1.0*StiffnessMatrix * U;
    if(lumping == 1)
        A = RHS ./ LumpedMass;
    else
        A = MassMatrix \ RHS;
    end
    t = t + delta_t;
    V = V + delta_t * A;
    %set boundary conditions (velocity)
    V(vDofsMoved) = Moved;
    %V(uDofsMoved) = Moved;
    U = U + delta_t * V;
    %set boundary conditions (displacement)
    U(uDofsFixed) = Fixed;
    U(vDofsFixed) = Fixed;
    iter = iter + 1;
    %StiffnessMatrix = CalcStiffnessMatrix(nn,Nodes+[U(1:nn) U(nn+1:2*nn)],Elements,C);
    %StiffnessMatrix(33,50);
end
%Calculate the stresses
disp([num2str(toc),' Calculate Stresses'])
Stresses = CalcStresses(Nodes,Elements,U,nn,C);
VonMisesStresses = zeros(size(Elements,1),1);
VonMisesStresses(:,1) = sqrt(Stresses(:,1).^2 + Stresses(:,2).^2 - Stresses(:,1).*Stresses(:,2) + 3*Stresses(:,3).^2);
%plot updated mesh
PlotMesh(Nodes+[U(1:nn) U(nn+1:2*nn)],Elements,c,L,VonMisesStresses,ax);
disp(['time: ',num2str(t)])
disp(['number of iterations: ', num2str(iter)])
disp(['norm of solution u: ', num2str(norm(U))])
disp(['maximum von mises stress: ', num2str(max(VonMisesStresses))])
disp([num2str(toc),' Calculation finished'])
end

