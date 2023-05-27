% Simulation zero dynamics of 3 link biped in the sagittal plane
%
% Inputs:
%       t: time
%       z: cyclic variables [q1, dq1]
%       a: bezier coefficient for q2 - a(1:5), and q3 - a(6:10)
%       s_params: [z_min, z_max] - max and min angles in gait
%       
% Outputs:       
%       dz = [dq1, ddq1]
%
function dz = func_zero_dynamics(t,z,a,s_params)

z_min = s_params(1); 
z_max = s_params(2); 
delq = z_max - z_min;

q1 = z(1); dq1 = z(2);

% Computes gait timing variable s
% 
% Inputs: 
%       q1: cyclic variable
%       q1_min: minimum value for q1 during gait
%       q1_max: max value for gait
%
% Output: 
%       s: gait timing variable
%
s = func_gait_timing(q1,z_min,z_max);

alpha2 = a(1:5);
alpha3 = a(6:10);

% find based on s and bezier definition
q2 = bezier(s,4,alpha2);
q3 = bezier(s,4,alpha3);
dq2 = d_ds_bezier(s,4,alpha2)*dq1/delq;
dq3 = d_ds_bezier(s,4,alpha3)*dq1/delq;

q = [q1, q2, q3];
dq = [dq1, dq2, dq3];

% Get model paramters
[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

% Compute D,C,G,B matrices
% Inputs:
%       q
%       dq
%       params
%
% Output: [D,C,G,B]
%
[D,C,G,~] = func_compute_D_C_G_B(q,dq,params);


%%%%%%%%%%%%%%%%%% Apply dynamics partitioning to compute the zero dynamics
%%%%%%%%%%%%%%%%%% equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E.g., 
D1 = D(1, 1);
D2 = [D(1, 2), D(1,3)];
H1 = C(1,:)*dq' + G(1);

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

% Compute 1st partial derivative of Bezier polynomial
% Inputs:
%       s: gait timing variable
%       M: order of Bezier polynomial, in this case 4th order
%       alpha: Coefficients for the polynomial, need M+1
%
% Outputs:  
%       db/ds

db_ds2 = d_ds_bezier(s,length(alpha2)-1,alpha2);
db_ds3 = d_ds_bezier(s,length(alpha3)-1,alpha3);

beta2 = [db_ds2; db_ds3]/delq;

%------------------------------------------------------------------------%

dz(1) = z(2);
dz(2) = -H1/D1;

dz = [dz(1), dz(2)]';
end