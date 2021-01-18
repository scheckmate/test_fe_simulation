function [ N, dNdxi ] = LagrangeBasis( dimension,pt )
%Calculate shape functions and gradient of shape functions for line and
%triangle respectively

if(dimension == 1)
    N = [1-pt; pt];
    dNdxi = [-1;1];
elseif(dimension == 2)
    N = [1-pt(1)-pt(2); pt(1); pt(2)];
    dNdxi = [-1 -1;
              1  0;
              0  1];
else
    N = 0;
    dNdxi = 0;
end
end
