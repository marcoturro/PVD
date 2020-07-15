function c = CfromI(DI)
% DI = |I - I0|
%%%%%%%%%---- For Live Stream ----%%%%%%%%%
% Goodness of fit:
%   SSE: 1.984e-32
% % Coefficients:
%        p1 =   -8.169e-06;
%        p2 =     0.004836;
%        p3 =   5.972e-17;
% c = p1*DI.^2 + p2*DI + p3;
%%%%%%%%%---- For Full Photos ----%%%%%%%%%
% Coefficients:
       p1 =  -1.408e-05;
       p2 =    0.005712;
       p3 =  -1.915e-16;
%Goodness of fit:
%  SSE: 6.749e-32
c = p1*DI.^2 + p2*DI + p3;
end