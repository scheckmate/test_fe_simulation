function [ C,E0,nu0,rho ] = GetMaterial( )
%Return some material dependent values
%Matrix C is for plane stress
%kg mm ms
%E0 = 210.0;
%nu0 = 0.3;
E0 = 70.0;
nu0 = 0.2;
rho = 7.81e-6;

C = E0/(1-nu0^2) *   [ 1   nu0    0;
                       nu0  1     0;
                       0    0 (1-nu0)/2 ];

end