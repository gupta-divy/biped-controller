
% Simulation of a single step of ZD using preimpact conditions applies 
%   impact first then simulates zero dynamics
%
% Inputs: f = [q10, dq10, alpha(3-5)_q2, alpha(3-5)_q3]
%       q10: pre-impact inital angle for q1
%       dq10: pre-impact inital velocity for dq1
%       alpha(3-5)_q2: 
%                   3rd to 5th Bezier coefficient for q2
%       alpha(3-5)_q3: 
%                   3rd to 5th Bezier coefficient for q3
%
% Outputs: 
%       t_sol - time (s) of zero dynamics
%       z_sol - [q1, dq1] post impact dynamics
%
function [t_sol, z_sol] = sim_zero_dynamics(f)

z_min = -f(1);  % To maintain symmetry min angle should be -ve of pre impact angle
z_max = f(1); % Pre impact angle should be max for gait

delq = z_max - z_min;

s_params = [z_min, z_max];

alpha2 = [-f(5)+2*f(3), -f(4)+2*f(3), f(3:5)];
alpha3 = [-f(8)+2*f(6), -f(7)+2*f(6), f(6:8)];

a = [alpha2, alpha3];

% Pre imapct conditions
z_minus = [f(1),f(2)];

%%%% Need all states to apply impact
%
% Mapping from z to x on boundary pg 140, 141

q2_minus = alpha2(5);            %q2- = alpha(M)
q3_minus = alpha3(5);            %q3- = gamma(M)

% d_dot- = M*(alpha(M) - alpha(M-1))*theta_dot-/(delta_theta)
dq2_minus = 4*(alpha2(5)-alpha2(4))*z_minus(2)/delq;
dq3_minus = 4*(alpha3(5)-alpha3(4))*z_minus(2)/delq;

x_minus = [z_minus(1), q2_minus, q3_minus, z_minus(2), dq2_minus, dq3_minus];

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
[x_plus, ~] = func_impact_map(x_minus);

% Intial condition for ODE45
z_plus = [x_plus(1), x_plus(4)];

%------------------------------------------------------------------------%

%%%% Simulation

tstart = 0; tfinal = 2;                    %max time per swing

options = odeset('Events',@events,'RelTol',1e-6,'AbsTol',1e-6);

% Integrate zero dynamics
% Inputs:
%       t: time
%       z: cyclic variables [q1, dq1]
%       a: bezier coefficient for q2 (1:5) and q3 (6:10)
%       s_param: [q1_min, q1_max]
%
% Outputs:
%       dz = [dq1, ddq1]
[t_sol, z_sol] = ode45(@(t,z) func_zero_dynamics(t,z,a,s_params), [tstart tfinal], z_plus, options);

%------------------------------------------------------------------------%

% Event function
% Inputs:
%       z: cyclic variables [q1, dq1]
%
    function [limit,isterminal,direction] = events(~,z)
        
        q1 = z(1);
        
        % Gait timing variable
        %
        % Inputs: [q1, q1_min, q1_max]
        %
        s = func_gait_timing(q1,z_min,z_max);
        
        % Saturate gait timing variable
        if s>= 1
            s = 1;
        elseif s<=0
            s = 0;
        end
        
        limit = s-1;
        isterminal = 1;    	% Halt integation
        direction = [];     %The zero can be approached from either direction
        
    end

end