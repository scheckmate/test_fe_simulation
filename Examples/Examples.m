function [ ] = Examples( )

%% symbolic integration in matlab
% symbolic integration is also possible in matlab but we are interested in
% numerical integration at the moment
disp('-----------------------------------')
disp('matlab symbolic')
syms t
disp(int(t))
disp(int(t,0,1))

%% numerical integration in matlab
disp('-----------------------------------')
disp('matlab numerical')
fun = @(x) x;
result1 = integral(fun,0,1);
disp(result1)
result2 = integral(fun,3,5);
disp(result2)

%% numerical integration using integration points and weigths from gaussian quadrature
disp('-----------------------------------')
disp('numerical integration using gaussian quadrature')
quad1 = NumIntegration(1,1);
fun = @(x) x;
result1 = 0;
for i=1:size(quad1,1)
    result1 = result1 + fun(quad1(i,1))*quad1(i,2);
end
disp(result1)

disp('numerical integration using gaussian quadrature and LagrangeBasis')
a = 3;
b = 5;
result1 = 0;
for i=1:size(quad1,1)
    point = quad1(i,1);
    [N,dNdxi] = LagrangeBasis(1,point);
    phi = [a b] * N;
    dphi = [a b] * dNdxi;
    result1 = result1 + fun(phi)*dphi*quad1(i,2);
end
disp(result1)

%% mass matrix
rho = 1;
elements = [1, 2; 2, 3];
nodes = [1.0; 3.0; 6.0];
numint = NumIntegration(1,2);
intpts = numint(:,1);
intweights = numint(:,2);
num_nodes = size(nodes,1);
massmatrix = zeros(num_nodes);
for e=1:size(elements,1)
    sctr = elements(e,:);
    for qtr=1:size(intpts,1)
        pt = intpts(qtr);
        weight = intweights(qtr);
        [N, dN] = LagrangeBasis(1, pt);
        dphi = nodes(sctr,:)' * dN;
        massqtr = rho * (N * N') * dphi * weight;
        massmatrix(sctr,sctr) = massmatrix(sctr,sctr) + massqtr;
    end
end

massmatrix
lumpedmass = sum(massmatrix,2)
% [L, U] = lu(massmatrix)
% L*U
%randRHS = rand(num_nodes,1)
% RHS = 0.5 * ones(num_nodes,1);
% x1 = massmatrix \ RHS
% x2 = RHS ./ lumpedmass

% massmatrix * x1
% lumpedmass .* x2
end

