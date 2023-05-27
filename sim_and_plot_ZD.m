
% Include util and autogen folders
set_path

%   f = [q10, dq10, alpha(3-5)_q2, alpha(3-5)_q3]
%       q10: pre-impact inital angle for q1
%       dq10: pre-impact inital velocity for dq1
%       alpha(3-5)_q2:
%                   3rd to 5th Bezier coefficient for q2
%       alpha(3-5)_q3:
%                   3rd to 5th Bezier coefficient for q3
%
% Manually tuned (working):
f = [   -0.3217   -1.2669   -0.0375    0.8621    0.6537    2.3562    3.2639    2.7826]; 
% (relatively good)
% f = [-0.2996 -2.2689 0.5571 1.2217 0.5992 2.9883 2.8947 2.9671];

z_tot = [];

n = 5;

for i = 1:n
    
    % Simulation of a single step of ZD to using preimpact conditions
    %   applies impact first then simulates zero dynamics
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
    [t_sol, z_sol] = sim_zero_dynamics(f);
    
    f(1:2) = z_sol(end,:);
    
    z_tot = [z_tot; z_sol];
    
end

figure
plot(z_tot(1,1),z_tot(1,2),'x'), hold on
plot(z_tot(:,1),z_tot(:,2))
hold off, grid on
title('Phase portrait of q_1 vd dq_1')
xlabel('q_1 (rads)')
ylabel('dq_1 (rads/s)')
legend('Start','Trajectory')