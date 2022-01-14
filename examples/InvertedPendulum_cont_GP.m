
%------------------------------------------%
% Inverted Pendulum 
  
% "Autogenerating microsecond solvers for nonlinear MPC: A tutorial
% using ACADO integrators", Quirynen, 2015

% typical configuration: 1) N=80, Ts=Ts_st=0.025, no shifting 2) N=40,
% Ts=Ts_st=0.05, shifting

%------------------------------------------%


%% Dimensions

nx=4;  % No. of differential states
nu=1;  % No. of controls
nz=0;  % No. of algebraic states
ny=5; % No. of outputs
nyN=4; % No. of outputs at the terminal point
npODE=5; % No. of model parameters
nc=0; % No. of general constraints
ncN=0; % No. of general constraints at the terminal point
nbx = 1; % No. of bounds on states
nbu = 1; % No. of bounds on controls

% GP
T = 200; %200 % No. of training points for the GP
xGP = 5; %5 % Dimension of training points for the GP
yGP = 2; %2 % # of target of the GP
npGP = T*xGP+yGP*(T+xGP); %X, ny*(alpha_i, l_i)

% state and control bounds
nbx_idx = 1; % indexs of states which are bounded
nbu_idx = 1; % indexs of controls which are bounded

% full params vector dimension
np = npODE+npGP;

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);   % differential states
controls = SX.sym('controls',nu,1); % control input
alg      = SX.sym('alg',nz,1);      % algebraic states
params   = SX.sym('paras',np,1); % parameters + GP params
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable

%% Dynamics

M = params(1);%1; 
m = params(2);%0.1;
l = params(3);%0.8; 

ode_flag = params(4);
gp_flag = params(5);

g = 9.81;

p=states(1);
theta=states(2);
v=states(3);
omega=states(4);  
u=controls(1);

a=-m*l*sin(theta)*omega^2+m*g*cos(theta)*sin(theta)+u;
b=-m*l*cos(theta)*sin(theta)*omega^2+u*cos(theta)+(M+m)*g*sin(theta);
c=M+m-m*(cos(theta))^2;

v_dot = a/c;
omega_dot = b/(l*c);

acc_fun = Function('acc_fun', {states,controls,params,alg},{[v_dot;omega_dot]});

%% GP definitions
addpath([pwd,'\gp_regression'])

% GP test point
x_star = [states',controls'];

% GP parametric functions
X = [];
for i = 1:xGP
    X = [X,params(npODE+(i-1)*T+1:npODE+i*T)];
end
X = X+SX.zeros(T,xGP);

idx_1 = npODE+xGP*T;
alpha_1 = params(idx_1+1:idx_1+T)+SX.zeros(T,1);
l_1 = params(idx_1+T+1:idx_1+T+xGP)'+SX.zeros(1,xGP);

idx_2 = idx_1+T+xGP;
alpha_2 = params(idx_2+1:idx_2+T)+SX.zeros(T,1);
l_2 = params(idx_2+T+1:idx_2+T+xGP)'+ SX.zeros(1,xGP);

% define GP for linear velocity
k_star_1 = RBF_cov_line_casadi(X, x_star, l_1);
psi_hat_1 = k_star_1*alpha_1;

% define GP for rotational velocity
k_star_2 = RBF_cov_line_casadi(X, x_star, l_2);
psi_hat_2 = k_star_2*alpha_2;

% define the GP values on the states
psi_bar = [0;...psi_hat_1*Ts_st/2; ... ...
          0;...psi_hat_2*Ts_st/2; ... ...
          psi_hat_1;...
          psi_hat_2...
         ];
% psi_bar = [psi_hat_1;...
%            psi_hat_2...
%          ];
       
settings.gp_generation = 'continuous';

%% explicit ODE RHS

x_dot=[ v;...
        omega;...
        v_dot*ode_flag;...
        omega_dot*ode_flag] + ...
        gp_flag*psi_bar;
 
% algebraic function
z_fun = [];                   

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;        
     
%% Objectives and constraints

% inner objectives
h = [p;theta;v;omega;u];
hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality constraints
general_con = [];
general_con_N = [];

%% NMPC discretizing time length [s]

Ts_st = 0.025; % shooting interval time
