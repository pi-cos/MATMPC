%% GP data save

num_points = 200;

% time
t = time';

% states
p       = state_sim(:,1)';
theta   = state_sim(:,2)';
v       = state_sim(:,3)';
omega   = state_sim(:,4)';

% controls
F = controls_MPC(:,1)';

if strcmp(settings.gp_generation,'discrete')

    % one-step-pred
    v_pred      = one_step_pred(:,3)';
    omega_pred  = one_step_pred(:,4)';

    % input: [p,theta,v,omega,V], [n x D]
    X = [p(1:end-1);theta(1:end-1);v(1:end-1);omega(1:end-1);F(1:end)]';

    % target [Delta v, Delta omega], [n x 2]
    if strcmp(opt.save_target,'disc_black_box')
        y = [v(2:end)-v(1:end-1);omega(2:end)-omega(1:end-1)]';
    elseif strcmp(opt.save_target,'disc_grey_box')
        y = [v(2:end)-v_pred;omega(2:end)-omega_pred;]';
    else
        error('Chosse correct saving options.')
    end

elseif strcmp(settings.gp_generation,'continuous')

%     % one-step-pred
%     v_dot      = one_step_pred(:,3)';
%     omega_dot  = one_step_pred(:,4)';

    % input: [p,theta,v,omega,V], [n x D]
%     X = [p(1:end-1);theta(1:end-1);v(1:end-1);omega(1:end-1);F(1:end)]';
    X = [p(1:end-1);theta(1:end-1);v(1:end-1);omega(1:end-1);F(1:end)]';

    % target [Delta v, Delta omega], [n x 2]
    if strcmp(opt.save_target,'cont_black_box')
        y = [acc_true(:,1),acc_true(:,2)];
    elseif strcmp(opt.save_target,'cont_grey_box')
        y = [acc_true(:,1)-acc_pred(:,1),acc_true(:,2)-acc_pred(:,2)];
    elseif strcmp(opt.save_target,'cont_grey_box_fromVel')
        acc_1_true_from_vel = diff(state_sim(:,3))/Ts;
        acc_2_true_from_vel = diff(state_sim(:,4))/Ts;
        y = [acc_1_true_from_vel-acc_pred(:,1),acc_2_true_from_vel-acc_pred(:,2)];
    elseif strcmp(opt.save_target,'cont_black_box_fromVel')
        acc_1_true_from_vel = diff(state_sim(:,3))/Ts;
        acc_2_true_from_vel = diff(state_sim(:,4))/Ts;
        y = [acc_1_true_from_vel,acc_2_true_from_vel];
    else
        error('Chosse correct saving options.')
    end

end

% concatenate in a table
T = array2table([t(1:end-1),X,y],...
    'VariableNames',{'t','x_1','x_2','x_3','x_4','x_5','y_1','y_2'});


% write table to csv
writetable(T(1:num_points,:),[pwd,'\gp_regression\GPR_PY\data\',name_csv,'.csv'])
    