function [ A ] = NumIntegration(dimension,order)
%Return integration points and weights for line (0 1) and
%triangle (0,0) - (1,0) - (0,1) respectively

% 1d segment, 1. column integration points, 2. column integration weights
% order 0 and 1
int_1d_1 = [0.5, 1.0];

% order 2 and 3
int_1d_3 = [(1.0-sqrt(1.0/3.0))/2.0  0.5;
    (1.0+sqrt(1.0/3.0))/2.0  0.5];

% order 4 and 5
int_1d_5 = [(1.0-sqrt(3.0/5.0))/2.0  5.0/18.0;
    0.5                     4.0/9.0;
    (1.0+sqrt(3.0/5.0))/2.0  5.0/18.0];

% order 6 and 7
int_1d_7 = [(1.0-sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)))/2.0  (18.0-sqrt(30.0))/72.0;
    (1.0-sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)))/2.0  (18.0+sqrt(30.0))/72.0;
    (1.0+sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)))/2.0  (18.0+sqrt(30.0))/72.0;
    (1.0+sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)))/2.0  (18.0-sqrt(30.0))/72.0];

% 2d triangle, 1. column integration points x, 2. column integration points
%y, 3. column integration weights
% order 0 and 1
int_2d_1 = [1.0/3.0, 1.0/3.0, 1.0/2.0];

% order 2
int_2d_2 = [1.0/6.0  1.0/6.0  1.0/6.0;
    4.0/6.0  1.0/6.0  1.0/6.0;
    1.0/6.0  4.0/6.0  1.0/6.0];
% order 3 and 4
int_2d_4 = [0.816847572980459  0.091576213509771  0.054975871827661;
    0.091576213509771  0.816847572980459  0.054975871827661;
    0.091576213509771  0.091576213509771  0.054975871827661;
    0.108103018168070  0.445948490915965  0.111690794839005;
    0.445948490915965  0.108103018168070  0.111690794839005;
    0.445948490915965  0.445948490915965  0.111690794839005];

if(dimension == 1)
    switch order
        case {0,1}
            A=int_1d_1;
        case {2,3}
            A=int_1d_3;
        case {4,5}
            A=int_1d_5;
        case {6,7}
            A=int_1d_7;
        otherwise
            A=0;
    end
elseif (dimension == 2)
    switch order
        case {0,1}
            A=int_2d_1;
        case 2
            A=int_2d_2;
        case {3,4}
            A=int_2d_4;
        otherwise
            A=0;
    end
else
    A=0;
end

end

