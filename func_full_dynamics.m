% Sagittal dynamics of 3 link biped
%
% Inputs:
%       t: time
%       x:
%           q1
%           q2
%           q3
%           dq1
%           dq2
%           dq3
%       alpha: Bezier coefficients for q2 and q3
%           alpha2 (1st to 5th coefficients)
%           alpha3 (1st to 5th coefficients)
%       s_params: for gait timing
%           q1_min
%           delq
%
% Outputs:
%       dx:
%           dq1
%           dq2
%           dq3
%           dd1
%           ddq2
%           ddq3
%
function dx = func_full_dynamics(t,x,alpha,s_params)

q = x(1:3);
dq = x(4:6);

[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

% Get D,C,G matrices
% Inputs:
%       q = [q1, q2, q3]
%       dq = [dq1, dq2, dq3]
%       params = [r,m,Mh,Mt,l,g]
[D,C,G,B] = func_compute_D_C_G_B(q,dq,params);
% Defining fx and gx
fx = [dq; D\(-C*dq - G)];
gx = [zeros(3,2); D\B];


%%%%%%%%%%%%%%%%%% Write your closed-loop feedback here %%%%%%%%%%%%%%%%%%%

% Get feedback
% Inputs:
%       x: states [q1,q2,q3,dq1,dq2,dq3]
%       alpha: Bezier coefficients for q2 and q3
%           alpha2 (1st to 5th coefficients)
%           alpha3 (1st to 5th coefficients)
%       s_params: for gait timing
%           q1_min
%           delq
u = func_feedback(x,alpha,s_params);
%u = [0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dx = fx + gx*u;

end