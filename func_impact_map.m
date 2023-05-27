% Obtain impact map from pre impact conditions - x_minus
% Inputs:
%       x_plus: states right before impact
%               [q1, q2, q3, dq1, dq2, dq3]
%
% Output:
%       x_plus: states right after impact
%               [q1, q2, q3, dq1, dq2, dq3]
%       F2: impact forces, x and y components
%

function [x_plus, F2] = func_impact_map(x_minus)

% Seperate states
q_minus = x_minus(1:3);
dq_minus = x_minus(4:6);
% Get params
[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

% Get matrices to compute impact
[De,E,dY_dq] = func_compute_De_E_dY_dq(q_minus, dq_minus, params);

dY_dq = zeros(2,3);

%Impact map from pg56 (3.20)

% Switching matrix that dictates how variable change roles post impact
R = [1, 1, 0;...
    0, -1, 0;...
    0, -1, 1]; 

delta_F2 = -(E*(De\E.'))\E*[eye(3); dY_dq];

% new angular angles
x_plus(1:3) = (R*x_minus(1:3)')';

% new angular velocities
x_plus(4:6) = [R, zeros(3,2)]*((De\E.')*delta_F2 + [eye(3); dY_dq])*x_minus(4:6).';

F2 = delta_F2;

end