% This .m file is used to symbolically derive:
%           the dynamics of the 3 link biped
%           impact map
%           controller
%           zero dynamics
% and write them onto function files
%
% Biped model: 
%   D, C, G matrices are found using a defined forward position kinematics,
%   B matrix must be defined manually 
%
% Impact map:
%   De, E. are derived using extended coordinates: p_e = [p_h; p_v]
%
% Controller:
%   The L2fh and LgLfh matrices used in feedback linearization are
%   symbolically derived
%
% Zero dynamics:
%   Vectors used in zero dynamics (eta2) are also derived
%

%-------------------------------------------------------------------------%
%%%% DCG matrices

syms q1 q2 q3 p_h p_v dp_h dp_v dq1 dq2 dq3 real
syms r m Mh Mt l g real

% Define parameters in a vector
params = [r,m,Mh,Mt,l,g];

% Include the util and autogen folder to use write_symbolic_term_to_mfile.m
% and export outputs to autogen folder
set_path

%Mh -  mass of hip, Mt - mass of torso, m - mass of legs
%l - length from hip to torso, r - length of legs

% Defining generalized coordinates:
% Angular positions:
%           q1: stance leg (absolute, w.r.t. y axis of 
%           q2: swing leg (relative to q1)
%           q3: torso (relative to q1)
% Angular velocities dq/dt:
%           dq1: stance leg 
%           dq2: swing leg 
%           dq3: torso  
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];

% q1 is cyclic, and negative pre-impact using convention provided in the
% figure

% Forward Kinematics - position of point masses
% hip
pMh = r*[sin(-q1); cos(-q1)];

% torso
pMt = pMh + l*[sin(-q1+pi-q3); cos(-q1+pi-q3)];

% stance leg
pm1 = r/2*[sin(-q1); cos(-q1)];

% swing leg
pm2 = pMh + r/2*[sin(q1+q2); -cos(q1+q2)];

% center of mass
pcm = simplify((m*pm1 + m*pm2 + Mh*pMh + Mt*pMt)/(m+m+Mh+Mt));

% end of swing leg
P2 = pMh + r*[sin(q1+q2); -cos(q1+q2)];


% Write positions to a file
% Inputs:
%       q
%       dq
%       params
%
% Outputs: Position vectors with x and y coordinates of position
%       pMh
%       pMt
%       pm1
%       pm2
%       pcm
%       P2
%
write_symbolic_term_to_mfile(q,dq,params,pMh,pMt,pm1,pm2,pcm,P2)

% Velocities - found by taking partial derivative w.r.t. q, then multiply
% by dq/dt
vMh = jacobian(pMh,q)*dq;
vMt = jacobian(pMt,q)*dq;
vm1 = jacobian(pm1,q)*dq;
vm2 = jacobian(pm2,q)*dq;
vcm = jacobian(pcm,q)*dq;

% Write velocities to a file
% Inputs:
%       q
%       dq
%       params
%
% Outputs: Velocity vectors with x and y coordinates of velocity
%       vMh
%       vMt
%       vm1
%       vm2
%       vcm
%       
write_symbolic_term_to_mfile(q,dq,params,vMh,vMt,vm1,vm2,vcm)


% Kinetic energy
%%%%%%%%%%%%%%%%%% Compute kinetic energy of each component here %%%%%%%%%%
K_Mh = 0.5*Mh*(vMh'*vMh);
K_Mt = 0.5*Mt*(vMt'*vMt);
K_m1 = 0.5*m*(vm1'*vm1);
K_m2 = 0.5*m*(vm2'*vm2);
% Total KE
K = K_m1+K_m2+K_Mt+K_Mh;

% Potential energy
%%%%%%%%%%%%%%%%%% Compute potnetial energy of each component here %%%%%%%%
V_Mh = Mh*g*pMh(2);
V_Mt = Mt*g*pMt(2);
V_m1 = m*g*pm1(2);
V_m2 = m*g*pm2(2);
% Total PE
V = V_m1 + V_Mh + V_Mt + V_m2;

% Inertia matrix
%%%%%%%%%%%%%%%%%% Compute D matrix here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = jacobian(jacobian(K,dq).',dq); D=simplify(D);

% Coriolis matrix
%%%%%%%%%%%%%%%%%% Coriolis matrix here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = max(size(q));
syms C
for k = 1:N
    for j = 1:N
        C(k,j) = 0*g;
        for i = 1:N
            C(k,j) = C(k,j)+0.5*(diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)))*dq(i);
        end
    end
end
C = simplify(C);

% Gravity matrix
%%%%%%%%%%%%%%%%%% Compute Gravity term here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = jacobian(V,q).';

%%%%%%%%%%%%%%%%%% What is control input matrix? %%%%%%%%%%%%%%%%%%%%%%%%%%
B = sym([0 0;1 0;0 1]);

% Write 3 link model to file
% Inputs:
%       q
%       dq
%       params
%
% Outputs: 
%       D: Inertia matrix
%       C: Coriolis matrix
%       G: Gravity matrix
%       B: 
%
write_symbolic_term_to_mfile(q,dq,params,D,C,G,B)

%-------------------------------------------------------------------------%
%%%% Impact map

% Using same psotion vectors as above, but taking partial with respect to qe
% instead 

% Extended configuration variables
p_e = [p_h; p_v];
qe = [q; p_h; p_v];
dqe = [dq; dp_h; dp_v];

% Extended position
pMh_e = pMh + p_e;
pMt_e = pMt + p_e;
pm1_e = pm1 + p_e;
pm2_e = pm2 + p_e;
P2e = P2 + p_e;

% Extended velocities
vMh_e = jacobian(pMh_e,qe)*dqe;
vMt_e = jacobian(pMt_e,qe)*dqe;
vm1_e = jacobian(pm1_e,qe)*dqe;
vm2_e = jacobian(pm2_e,qe)*dqe;

K_Mhe = 0.5*Mh*(vMh_e'*vMh_e); 
K_Mte = 0.5*Mt*(vMt_e'*vMt_e);
K_m1e = 0.5*m*(vm1_e'*vm1_e);
K_m2e = 0.5*m*(vm2_e'*vm2_e);

Ke = K_m1e + K_Mhe + K_Mte + K_m2e;

% Extended inertia matrix
De = jacobian(jacobian(Ke,dqe).',dqe); De=simplify(De);

E = jacobian(P2e,qe);

% Partial of any point on biped, hip chosen in this case
dY_dq = jacobian(pMh_e,q);

% Write impact map to a file
% Inputs:
%       q
%       dq
%       params
%
% Outputs: Matrices needed to compile impact map
%       De: Extended inertia matrix
%       E:
%       dY_dq:
%       
write_symbolic_term_to_mfile(q,dq,params,De,E,dY_dq)


%-------------------------------------------------------------------------%
%%%% For controller
% Vector fields
fx = [dq; D\(-C*dq - G)];
gx = [sym(zeros(3,2));D\B];

% Bezier poly - needed for output function
syms s delq
%s = (q1 - q1_plus)/delq; 
%delq = q1_minus - q1_plus 
%ds/dt = dq1/delq; ds/dq1 = 1/delq;

syms a21 a22 a23 a24 a25 
syms a31 a32 a33 a34 a35

a2 = [a21 a22 a23 a24 a25];
a3 = [a31 a32 a33 a34 a35];
M = 4;

b2 = 0; b3 = 0;
for k = 0:M
    b2 = b2 + a2(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

for k = 0:M
    b3 = b3 + a3(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

% Defining outputs

h = [q2 - b2; q3 - b3];

% y_dot = Lfh = dh/dx*fx - independent of gx*u since relative degree is 2
% However, h is a function of (s,q2,q3), not q1 directly, so the following
% is used:
% dh/dq1 = dh/ds*ds/dq1 = dh/ds*1/delq
%
% Temporary variable that multiples the 1st column with 1/delq
temp = sym(eye(6)); temp(1)  = 1/delq;

Lfh = jacobian(h,[s; q2; q3; dq])*temp*fx;

dLfh = jacobian(Lfh,[s; q2; q3; dq])*temp;

% Write matrix used in feedback linearization - d/dx(Lfh) to file
% Inputs:
%       s = (q1 - q1_plus)/delq: gait timing variable
%       delq = q1_minus - q1_plus: difference in cyclic variable during gait 
%       dq1
%       params: 
%       a2: bezier coefficents (1st - 5th) for q2
%       a3: bezier coefficents (1st - 5th) for q3
%
% Outputs:
%       dLfh: partial of Lfh, to be used to compute L2fh and LgLfh
%
write_symbolic_term_to_mfile([s,delq],dq1,[a2,a3],dLfh);

%-------------------------------------------------------------------------%
%%%% For Zero Dynamics

db_ds2 = 0;
for k = 0:M-1
    db_ds2 = db_ds2 + (a2(1,k+2)-a2(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

db_ds3 = 0;
for k = 0:M-1
    db_ds3 = db_ds3 + (a3(1,k+2)-a3(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

partial_db_ds2 = jacobian(db_ds2,s)*dq1/delq;

partial_db_ds3 = jacobian(db_ds3,s)*dq1/delq;

beta1 = [partial_db_ds2; partial_db_ds3]*dq1/delq;

eta2 = jacobian(K,dq1);

write_symbolic_term_to_mfile(s,[dq1, delq],[a2, a3],beta1)

write_symbolic_term_to_mfile(q,dq,params,eta2)
