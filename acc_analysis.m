%%

restoredefaultpath; clear all; %close all; clc;

%%

addpath('model_src')
addpath('mex_core')

%%

load sim_res.mat

Ts_orig = settings.Ts_st;
Ts_factor = 1;
Ts = Ts_orig/Ts_factor;
Tf = 10;

[input, data] = InitData(settings);

mem.num_steps = 2;
mem.h=Ts/mem.num_steps;

%% computations

acc_true = [];
acc_nom = [];

vel_true = [];
vel_nom = [];

int_vel_erk = [];
int_vel_eul = [];

state_sim = input.x0';

for ii = 1:round(Tf/Ts)-Ts_factor
    % Simulate REAL system dynamics
    sim_input.x = state_sim(end,:).';
    sim_input.u = controls_MPC(round(ii/Ts_factor));
    sim_input.z = input.z(:,1);
    sim_input.p = input.od(:,1);% -> real is [1;0.1;0.5]
    sim_input.p(1:settings.npODE) = [1;0.1;0.5;1;0];

    acc_true_curr = acc_fun('acc_fun',sim_input.x,sim_input.u,sim_input.p,sim_input.z );
    acc_true_curr_sym = acc(sim_input.x,sim_input.u,sim_input.p(1:3));
    acc_true = [acc_true;acc_true_curr',acc_true_curr_sym(3:4)'];

    int_vel_erk = [int_vel_erk;...
                   ERK(sim_input.x,sim_input.u,sim_input.p(1:3),Ts,'acc')'];

    int_vel_eul = [int_vel_eul;...
                   Euler(sim_input.x,sim_input.u,sim_input.p(1:3),Ts,'acc')'];

    [xf, ~] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
    xf = full(xf);
    state_sim = [state_sim; xf'];
    vel_true = [vel_true; xf(3:4)'];


    % Simulate NOMINAL system dynamics
    sim_input.p = input.od(:,1);% -> nominal is [1;0.1;0.8]
    sim_input.p(1:settings.npODE) = [1;0.1;0.8;1;0];

    acc_nom_curr = acc_fun('acc_fun',sim_input.x,sim_input.u,sim_input.p,sim_input.z );
    acc_nom_curr_sym = acc(sim_input.x,sim_input.u,sim_input.p(1:3));
    acc_nom = [acc_nom;acc_nom_curr',acc_nom_curr_sym(3:4)'];

    [xf, zf] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
    xf = full(xf);
    vel_nom = [vel_nom; xf(3:4)'];
end

%%

figure
% subplot(211)
% hold on
% plot(acc_nom(2:end,1),'-')
% plot(diff(vel_nom(:,1))/Ts,'-.')
% title('linear acc')
% subplot(212)
hold on
plot(acc_true(2:end,2),'-')
plot(diff(vel_true(:,2))/Ts,'-.')
title('rotational acc')
legend('true acc','true vel diff')


%%

figure
hold on
plot(cumsum(acc_true(:,2))*Ts,'-')
plot(vel_true(:,2),':')
plot(int_vel_erk(:,4),'--')
plot(int_vel_eul(:,4),'-.')
title('rotational vel')
legend('cumsum','sym sys','erk','eul')

%%

diff_2 = get_first_derivative(int_vel_erk(:,4),Ts);
[q_3,q_dot_3,diff_3] = get_derivatives(int_vel_erk(:,2),Ts);

figure
hold on
plot(acc_true(4:end,2),'-')
plot(diff(int_vel_erk(3:end,4))/Ts,'--')
plot(diff_2,'--')
plot(diff_3,'--')
% plot(diff(int_vel_eul(3:end,4))/Ts,'-.') 
title('rotational acc')
legend('true acc','erk diff','erk diff 2','erk diff 3')%,'eul diff'

%%
% Euler(sim_input.x,sim_input.u,sim_input.p(1:3),Ts,'acc')
% ERK(sim_input.x,sim_input.u,sim_input.p(1:3),Ts,'acc')

%% integrators

% ERK(sim_input.x,sim_input.u,sim_input.p(1:3),Ts,'acc')

function X_new = ERK(X,U,P,Ts,f)
s  = 1; % No. of integration steps per sample interval
DT = Ts/s;
for j=1:s
       k1 = feval(f,X, U, P);
       k2 = feval(f,X + DT/2 * k1, U, P);
       k3 = feval(f,X + DT/2 * k2, U, P);
       k4 = feval(f,X + DT * k3, U, P);
       X_new=X+DT/6*(k1 +2*k2 +2*k3 +k4);
end
end

% Euler(sim_input.x,sim_input.u,sim_input.p(1:3),Ts,'acc')
function X_new = Euler(X,U,P,Ts,f)
X_new = X+Ts*feval(f,X, U, P);
end

%% acceleration

function x_dot = acc(states,controls,params)

p=states(1);
theta=states(2);
v=states(3);
omega=states(4);  

u=controls(1);

g = 9.81;

M = params(1);%1; 
m = params(2);%0.1;
l = params(3);%0.8; 

a=-m*l*sin(theta)*omega^2+m*g*cos(theta)*sin(theta)+u;
b=-m*l*cos(theta)*sin(theta)*omega^2+u*cos(theta)+(M+m)*g*sin(theta);
c=M+m-m*(cos(theta))^2;

v_dot = a/c;
omega_dot = b/(l*c);

x_dot=[ v;...
        omega;...
        v_dot;...
        omega_dot];
end

%% derivatives

function q_dot = get_first_derivative(q,Ts)

% fs = 1/Ts;
% fc = 10;
% [b,a] = butter(2,fc/(2*fs));
% q = filtfilt(b,a,q);

q_dot = (q(3:end)-q(1:end-2))/(2*Ts);

end

function [q,q_dot, q_ddot] = get_derivatives(q,Ts)

% Ts = time(2)-time(1);

% fs = 1/Ts;
% fc = 10;
% [b,a] = butter(2,fc/(2*fs));
% q = filtfilt(b,a,q);

q_dot = (q(3:end)-q(1:end-2))/(2*Ts);
q_ddot = (q_dot(3:end)-q_dot(1:end-2))/(2*Ts);
q = q(3:end-2);
q_dot = q_dot(2:end-1);
% time_diff = time(3:end-2);

% q_dot = diff(q)/Ts;
% q_ddot = diff(q_dot)/Ts;
% 
% q = q(3:end);
% q_dot = q_dot(2:end);
% 
% time_diff = time(3:end);

end