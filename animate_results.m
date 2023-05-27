%
% Input:
%       t: time(s) for the corresponding gait in x
%       x: [q1, q2, q3, dq1, dq2, dq3]
%
% angles are in radians, velocities are in rad/s
%

function animate_results(t,x)

% Time is always 0 at the begining of an ODE45 solution, so it is the most
% consistent way to seperate phases and steps
%
% Find index when time = 0:
ind0 = find(t == 0);
% Number of steps taken:
n = length(ind0);

% Add extra element to index to help with changing origin after each step
ind0(length(ind0)+1) = length(t)+1;

% Get model paramters
[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

%%% Estimate x axis limits for animation window
% Find step length at begining of gait
[~,~,~,~,~, P2_i] = func_compute_pMh_pMt_pm1_pm2_pcm_P2(x(1,1:3),x(1,4:6),params);
delta = 0.5;
xlim_min = P2_i(1) - delta;
xlim_max = -n*P2_i(1) + delta;

% Defining figure properties
fh = figure('Name','3 link biped model in the sagittal plane',...
    'Renderer','opengl',...
    'GraphicsSmoothing','on');
ah = axes('Box','on',...
    'XGrid','off',...
    'YGrid','off',...
    'DataAspectRatio',[1,1,1],...
    'PlotBoxAspectRatio',[1,1,1],...
    'Parent',fh);
xlabel(ah,'[m]');
ylabel(ah,'[m]');
xlim(ah,[xlim_min, xlim_max]);
ylim(ah,[-0.5 2]);

hold(ah,'off'); % To animate

for i = 1:length(t)
    
    q = x(i,1:3);
    dq = x(i,4:6);
    
    [pMh,pMt,pm1,pm2,~, P2] = func_compute_pMh_pMt_pm1_pm2_pcm_P2(q,dq,params);
    
    % If it is the 1st step, origin should be at (0,0)
    if i < ind0(2)
        L_step = [0; 0];
    % If step number > 1, origin must be shifted by step length of previous
    % step
    elseif sum(i == ind0) && i ~= 1 % True only if i matches one of the entries of ind0 except 1
        % Holds value until next index matches
        j = i-1;  % Index to find step length 
        [~,~,~,~,~, P2] = func_compute_pMh_pMt_pm1_pm2_pcm_P2(x(j,1:3),x(j,4:6),params);
        L_step = L_step + P2;  % Add to previous step length
    end
    
    % Clear figure 
    cla(ah);
    % Add dots to begining of stance leg and end of swing leg
    line(ah, [L_step(1),L_step(1)], [L_step(2),L_step(2)],'Color','r','Marker','.')
    line(ah, L_step(1)+[P2(1),P2(1)], L_step(2)+[P2(2),P2(2)],'Color','g','Marker','.')
    
    % Plot links as dashed lines
    line(ah, L_step(1)+[0,pMh(1)], L_step(2)+[0,pMh(2)],'Color','r','LineStyle','--')
    line(ah, L_step(1)+[pMh(1),pMt(1)], L_step(2)+[pMh(2),pMt(2)],'Color','k','LineStyle','--')
    line(ah, L_step(1)+[pMh(1),P2(1)], L_step(2)+[pMh(2),P2(2)],'Color','g','LineStyle','--')
    
    % Plot point mass as circles
    line(ah, L_step(1)+[pm1(1),pm1(1)], L_step(2)+[pm1(2),pm1(2)],'Color','r','Marker','o')
    line(ah, L_step(1)+[pm2(1),pm2(1)], L_step(2)+[pm2(2),pm2(2)],'Color','g','Marker','o')
    line(ah, L_step(1)+[pMh(1),pMh(1)], L_step(2)+[pMh(2),pMh(2)],'Color','k','Marker','o')
    line(ah, L_step(1)+[pMt(1),pMt(1)], L_step(2)+[pMt(2),pMt(2)],'Color','k','Marker','o')
    
    %legend([p1, p2, p3],{'stance leg','swing leg','torso'});

    drawnow limitrate
    %pause(0.02);
end