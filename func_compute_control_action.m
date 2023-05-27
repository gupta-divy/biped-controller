% Post compute control action
%
% Inputs: [z, a, s_params]
%       z: [q1, dq1]
%       a: [alpha2, alpha3] - 1st to 5th Bezier coefficient for q2 and q3
%       s_params: [q1_min, q1_max]
%
% Output:
%       u: control action
%
function u = func_compute_control_action(z,a,s_params)

% Seperate inputs
q1_min = s_params(1);
q1_max = s_params(2);
delq = q1_max - q1_min;

% Maps zero dynamics coordinates to full states using bezier polynomials
%
% Inputs:
%       z: [z, dz] = [q1, dq1]
%           stance leg angle and velocity
%       a: [alpha_2, alpha_3], bezier coefficients for q2 and q3
%           each has 5 coefficients for 4th order Bezier Poly
%       s_params: [q1_min, q1_max]
%           q1_min: minimum angle during gait
%           q1_max: maxmum angle during gait
%
% Outputs:
%       x: full dynamics [q, dq]
x = func_map_z_x(z,a,s_params);

q = x(1:3);
dq = x(4:6);

[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

%------------------------------------------------------------------------%

%%%% Partitioned dynamics

% Compute D,C,G,B matrices
% Inputs:
%       q
%       dq
%       params
%
% Output: [D,C,G,B]
%
[D,C,G,B] = func_compute_D_C_G_B(q,dq,params);      

% Bezier coefficients for q2;
alpha2 = a(1:5);
% Bezier coefficients for q3
alpha3 = a(6:10);

q1 = z(1);
dq1 = z(2);

D1 = D(1, 1);
D2 = D(1, 2:3);
D3 = D(2:3, 1);
D4 = D(2:3, 2:3);

H = C*dq'+ G;
H1 = H(1);
H2 = H(2:3);
% Computes gait timing variable s
% 
% Inputs: 
%       q1: cyclic variable
%       q1_min: minimum value for q1 during gait
%       q1_max: max value for gait
%
s = func_gait_timing(q1,q1_min,q1_max);   %normalized general coordinate

%------------------------------------------------------------------------%

% Get beta1 term to compute dz
%
% Inputs: s, [dq1, delq], [alpha2, alpha3]
%       s: gait timing variable
%       dq1: velocity of cyclic variable = z2
%       delq: difference between z_max and z_min in s, needed when
%               computing db/ds
%       alpha2: Bezier coefficient (1st to 5th) for q2
%       alpha3: Bezier coefficient (1st to 5th) for q3
%
% Output:
%       beta1: one of 2 matrices obtained when computing ddq_b (acceleration
%               of body coordinates) 
%       	y = h(x) = q_b - b(s), when y = 0,  q_b = b(s)
%           taking two time derivatives we get:
%       	ddq_b = d/ds(db(s)/ds*ds/dt)*ds/dt + db(s)/ds * d^2(s)/dt^2
%                 = beta1 + beta2*ddq_1
%           Note: d^2(s)/dt^2 = 1/delq*d^2(q1)/dt^2
%
%       	beta1 = d/ds(db(s)/ds*ds/dt)*ds/dt
%        	beta2 = db(s)/ds*1/delq
%
beta1 = func_compute_beta1(s, [dq1, delq], [alpha2, alpha3]);

db_ds2 = d_ds_bezier(s,4,alpha2);
db_ds3 = d_ds_bezier(s,4,alpha3);

beta2 = [db_ds2; db_ds3]/delq;

%------------------------------------------------------------------------%
ddq1 = (D1 + D2*beta2)\(-D2*beta1 - H1);
u = B(2:3,1:2)\((D3 + D4*beta2)*ddq1 + (D4*beta1 + H2));
end