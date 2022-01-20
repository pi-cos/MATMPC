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

state_sim = input.x0';

for ii = 1:round(Tf/Ts)-Ts_factor
    % Simulate REAL system dynamics
    sim_input.x = state_sim(end,:).';
    sim_input.u = controls_MPC(round(ii/Ts_factor)+1);
    sim_input.z = input.z(:,1);
    sim_input.p = input.od(:,1);% -> real is [1;0.1;0.5]
    sim_input.p(1:settings.npODE) = [1;0.1;0.5;1;0];

    acc_true_curr = acc_fun('acc_fun',sim_input.x,sim_input.u,sim_input.p,sim_input.z );
    acc_true = [acc_true;acc_true_curr'];

    [xf, ~] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
    xf = full(xf);
    state_sim = [state_sim; xf'];
    vel_true = [vel_true; xf(3:4)'];

    % Simulate NOMINAL system dynamics
    sim_input.p = input.od(:,1);% -> nominal is [1;0.1;0.8]
    sim_input.p(1:settings.npODE) = [1;0.1;0.8;1;0];

    acc_nom_curr = acc_fun('acc_fun',sim_input.x,sim_input.u,sim_input.p,sim_input.z );
    acc_nom = [acc_nom;acc_nom_curr'];

    [xf, zf] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
    xf = full(xf);
    vel_nom = [vel_nom; xf(3:4)'];
end

%%

figure
subplot(211)
hold on
plot(acc_nom(2:end,1),'-')
plot(diff(vel_nom(:,1))/Ts,'-.')
title('linear acc')
subplot(212)
hold on
plot(acc_nom(2:end,2),'-')
plot(diff(vel_nom(:,2))/Ts,'-.')
title('rotational acc')
legend('nominal acc','nominal vel diff')

% figure(2)
% subplot(211)
% hold on
% plot(cumsum(acc_nom(1:end,1))*Ts,'-')
% plot(vel_nom(:,1),'-.')
% title('linear acc')
% subplot(212)
% hold on
% plot(cumsum(acc_nom(1:end,2))*Ts,'-')
% plot(vel_nom(:,2),'-.')
% title('rotational acc')
% legend('nom acc int','nom vel')
