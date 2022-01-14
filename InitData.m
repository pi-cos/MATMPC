%% Initialize Data
function [input, data] = InitData(settings)

    nx = settings.nx;       % No. of differential states
    nu = settings.nu;       % No. of controls
    nz = settings.nz;       % No. of algebraic states
    ny = settings.ny;        % No. of outputs (references)    
    nyN= settings.nyN;       % No. of outputs at terminal stage 
    np = settings.np;        % No. of parameters (on-line data)
    nc = settings.nc;        % No. of constraints
    ncN = settings.ncN;      % No. of constraints at terminal stage
    N  = settings.N;         % No. of shooting points
    nbx = settings.nbx;      % No. of state bounds
    nbu = settings.nbu;      % No. of control bounds
    nbu_idx = settings.nbu_idx;  % Index of control bounds

    switch settings.model
                      
        case 'InvertedPendulum'
            input.x0 = [0;pi;0;0];    
            input.u0 = zeros(nu,1); 
            input.z0 = zeros(nz,1);
            ode_flag = 1;
            para0 = [1;0.1;0.8;ode_flag];  

            Q=repmat([10 10 0.1 0.1 0.01]',1,N);
            QN=[10 10 0.1 0.1]';

            % upper and lower bounds for states (=nbx)
            lb_x = -2;
            ub_x = 2;

            % upper and lower bounds for controls (=nbu)           
            lb_u = -50;
            ub_u = +50;
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];
            
        case 'InvertedPendulum_disc_GP'
            
            input.x0 = [0;pi;0;0];    
            input.u0 = zeros(nu,1); 
            input.z0 = zeros(nz,1);
            if strcmp(settings.mpc_model,'white_box_corr')
                l_p = 0.5;
            else
                l_p = 0.8;
            end
            if strcmp(settings.mpc_model,'black_box')
                ode_flag = 0;
            else
                ode_flag = 1;
            end
            para0_ode = [1;0.1;l_p;ode_flag];  

            Q=repmat([10 10 0.1 0.1 0.01]',1,N);
            QN=[10 10 0.1 0.1]';

            % upper and lower bounds for states (=nbx)
            lb_x = -2;
            ub_x = 2;

            % upper and lower bounds for controls (=nbu)           
            lb_u = -50;
            ub_u = +50;
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];
            
            % GP params
            if strcmp(settings.mpc_model,'grey_box') || strcmp(settings.mpc_model,'black_box')
                if strcmp(settings.mpc_model,'black_box')
                    gp_res_path = [pwd,'\gp_regression\GPR_PY\results_GP_ID\blackbox'];
                else %strcmp(settings.mpc_model,'grey_box')
                    gp_res_path = [pwd,'\gp_regression\GPR_PY\results_GP_ID\greybox'];
                end
                name = '';
                X = readmatrix([gp_res_path,'\X',name,'.csv']);
                alpha_1 = readmatrix([gp_res_path,'\alpha_1',name,'.csv']);
                l_1 = readmatrix([gp_res_path,'\l_1',name,'.csv'])';
                alpha_2 = readmatrix([gp_res_path,'\alpha_2',name,'.csv']);
                l_2 = readmatrix([gp_res_path,'\l_2',name,'.csv'])';
                
                X_line = [];
                for i =1:settings.xGP
                    X_line = [X_line;X(:,i)];
                end
                para0_GP = [X_line;alpha_1;l_1';alpha_2;l_2'];
            else
                para0_GP = zeros(settings.npGP,1);%[X_line;alpha_1;l_1';alpha_2;l_2'];
            end
            
%             para0_GP = zeros(settings.npGP,1);
            
            para0 = [para0_ode;para0_GP];

        case 'InvertedPendulum_cont_GP'
            input.x0 = [0;pi;0;0];
            input.u0 = zeros(nu,1);
            input.z0 = zeros(nz,1);
            ode_flag = 1;
            gp_flag = 1;
            para0_ode = [1;0.1;0.8;ode_flag;gp_flag];

            Q=repmat([10 10 0.1 0.1 0.01]',1,N);
            QN=[10 10 0.1 0.1]';

            % upper and lower bounds for states (=nbx)
            lb_x = -2;
            ub_x = 2;

            % upper and lower bounds for controls (=nbu)
            lb_u = -50;
            ub_u = +50;

            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];
            lb_gN = [];
            ub_gN = [];

            % GP params
            if gp_flag
%                 if strcmp(settings.mpc_model,'black_box')
%                     gp_res_path = [pwd,'\gp_regression\GPR_PY\results_GP_ID\blackbox'];
%                 else %strcmp(settings.mpc_model,'grey_box')
                    gp_res_path = [pwd,'\gp_regression\GPR_PY\results_GP_ID\cont_greybox'];
%                 end
                name = '';
                X = readmatrix([gp_res_path,'\X',name,'.csv']);
                alpha_1 = readmatrix([gp_res_path,'\alpha_1',name,'.csv']);
                l_1 = readmatrix([gp_res_path,'\l_1',name,'.csv'])';
                alpha_2 = readmatrix([gp_res_path,'\alpha_2',name,'.csv']);
                l_2 = readmatrix([gp_res_path,'\l_2',name,'.csv'])';
                
                X_line = [];
                for i =1:settings.xGP
                    X_line = [X_line;X(:,i)];
                end
                para0_GP = [X_line;alpha_1;l_1';alpha_2;l_2'];
            else
                para0_GP = zeros(settings.npGP,1);%[X_line;alpha_1;l_1';alpha_2;l_2'];
            end

            para0 = [para0_ode;para0_GP];

    end

    % prepare the data
    
    input.lb = repmat(lb_g,N,1);
    input.ub = repmat(ub_g,N,1);
    input.lb = [input.lb;lb_gN];
    input.ub = [input.ub;ub_gN];
            
    lbu = -inf(nu,1);
    ubu = inf(nu,1);
    for i=1:nbu
        lbu(nbu_idx(i)) = lb_u(i);
        ubu(nbu_idx(i)) = ub_u(i);
    end
                
    input.lbu = repmat(lbu,1,N);
    input.ubu = repmat(ubu,1,N);
    
    input.lbx = repmat(lb_x,1,N);
    input.ubx = repmat(ub_x,1,N);
        
    x = repmat(input.x0,1,N+1);  % initialize all shooting points with the same initial state
    u = repmat(input.u0,1,N);    % initialize all controls with the same initial control
    z = repmat(input.z0,1,N);    % initialize all algebraic state with the same initial condition
    para = repmat(para0,1,N+1);  % initialize all parameters with the same initial para
         
    input.x=x;           % (nx by N+1)
    input.u=u;           % (nu by N)
    input.z=z;           % (nz by N)
    input.od=para;       % (np by N+1)
    input.W=Q;           % (ny by N)
    input.WN=QN;         % (nyN by 1)
     
    input.lambda=zeros(nx,N+1);   % langrangian multiplier w.r.t. equality constraints
    input.mu=zeros(N*nc+ncN,1);   % langrangian multipliers w.r.t. general inequality constraints
    input.mu_u = zeros(N*nu,1);   % langrangian multipliers w.r.t. input bounds
    input.mu_x = zeros(N*nbx,1);  % langrangian multipliers w.r.t. state bounds
    
    %% Reference generation

    switch settings.model

        case {'InvertedPendulum','InvertedPendulum_disc_GP','InvertedPendulum_cont_GP'}

            data.REF=zeros(1,nx+nu);

                     
    end
    
end