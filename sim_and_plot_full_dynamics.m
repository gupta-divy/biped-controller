% Simulates full dynamics for a 3 link biped in the sagittal plane
% Needed "inputs" from optimizer:
%   f = [q10, dq10, alpha(3-5)_q2, alpha(3-5)_q3]
%       q10: pre-impact inital angle for q1
%       dq10: pre-impact inital velocity for dq1
%       alpha(3-5)_q2: 
%                   3rd to 5th Bezier coefficient for q2
%       alpha(3-5)_q3: 
%                   3rd to 5th Bezier coefficient for q3
%

%------------------------------------------------------------------------%

% Include util and autogen folders
set_path

% Ouput of optimizer
% The f vector below is actually walking down a slope - optimizer wasn't
% enforcing required constraints
% f = [-0.4821, -3.597, 0.2048, 0.1918, 0.3672, 0.6869, 2.3175, 2.2645, 2.8316];
% x_max = 0.2048; 


%   f = [q10, dq10, alpha(3-5)_q2, alpha(3-5)_q3]
%       q10: pre-impact inital angle for q1
%       dq10: pre-impact inital velocity for dq1
%       alpha(3-5)_q2: 
%                   3rd to 5th Bezier coefficient for q2
%       alpha(3-5)_q3: 
%                   3rd to 5th Bezier coefficient for q3
%

f = [   -0.3217   -1.2669   -0.0375    0.8621    0.6537    2.3562 3.2639    2.7826];
% f= [-0.3478,-3.0765,0.1520,0.3604,0.5711,2.8535,2.8806,3.2398];
%------------------------------------------------------------------------%

% extract bezier coeffients
alpha2 = [-f(5), -f(4), f(3:5)];
alpha3 = [-f(8)+2*f(6), -f(7)+2*f(6),f(6:8)];

a = [alpha2,alpha3];

% extract preimact states
q1_minus = f(1);
dq1_minus = f(2);

% maximum and minimum angles for q1 during a sigle gait
x_max = q1_minus;  
x_min = -q1_minus;      % negative to ensure symmetry within gait
delq = x_max - x_min;

%------------------------------------------------------------------------%
%%%% Impact map
% Need full states to get impact map

%Mapping from z to x on boundary pg 140, 141
q2_minus = alpha2(5);
q3_minus = alpha3(5);

%d_dot- = M*(alpha(M) - alpha(M-1))*theta_dot-/(delta_theta)
dq2_minus = 4*(alpha2(5)-alpha2(4))*dq1_minus/delq;
dq3_minus = 4*(alpha3(5)-alpha3(4))*dq1_minus/delq;

%I.C. [thetas; velocities]
x_minus = [q1_minus, q2_minus, q3_minus, dq1_minus, dq2_minus, dq3_minus];

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
[x_plus0, ~] = func_impact_map(x_minus);

s_params = [x_plus0(end,1), x_max];

%%%% Simulate
[t_tot,x_tot] = sim(x_plus0, a, s_params);

%------------------------------------------------------------------------%

% Plots and animation

plot_trajectories(t_tot,x_tot)

animate_results(t_tot,x_tot)

%------------------------------------------------------------------------%

% Simulater, wrapped in a function so event function can accept global
% variable
%
% Inputs:
%   x_plus0: [q0, dq0]
%           Intial conditions for ODE45
%       a: [alpha2, alpha3]
%           Bezier coefficients for q2 and q3
%       s_params: [q1_min, q1_max]      
%
% Outputs:
%   t_tot: times (in s) for all steps appended together, time = 0 at the 
%       begining of each step. This will help to sort steps in plots 
%   x_tot: [q, dq]
%       output of ODE45 for all steps appended together
%
function [t_tot,x_tot] = sim(x_plus0, a, s_params)

q1_min = s_params(1);
q1_max = s_params(2);

ti = 0; tf = 5;

options = odeset('Event',@event,'AbsTol',1e-6,'RelTol',1e-6);

n = 10;     % Number of desired steps the biped should take

% Define variables where time and states of solution will be appended
t_tot = []; x_tot = [];

for i = 1:n
    
    [t,x_sol] = ode45(@(t,x) func_full_dynamics(t,x,a,s_params), [ti, tf], x_plus0, options);
    
    % Apply impact again - only useful is n > 1
    [x_plus, ~] = func_impact_map(x_sol(end,:));
    
    % New initial condition
    x_plus0 = x_plus(end,:);
    
    if isempty(t_tot)
        t_end = 0;
    end
    
    % append
    t_tot = [t_tot; t_end + t];
    x_tot = [x_tot; x_sol];
    
    disp(['Step#...',num2str(i)]);
    
end

%------------------------------------------------------------------------%

% Event function - detect when impact happens
%
% Inputs: [t, x]
%
% Note: currently using gait timing variable instead of end of swing feet
%
    function [limits,terminal,direction] = event(~,x)
        
        q = x(1:3);
        dq = x(4:6);
        
        s = func_gait_timing(q(1),q1_min,q1_max);
        
        [r,m,Mh,Mt,l,g] = func_model_params;
        params = [r,m,Mh,Mt,l,g];
        
        [~,~, pm1, pm2, P2] = func_compute_pMh_pMt_pm1_pm2_pcm_P2(q,dq,params);
        [~,~,~, ~, vcm] = func_compute_vMh_vMt_vm1_vm2_vcm(q,dq,params);
        
        %limits = P2(2) <= 0.01 && pm2(1) > pm1(1) && vcm(2) < 0;
        
        if s>=1
            s=1;
        end
        
        limits = s-1;
        terminal = 1;
        direction = [];
    end

end