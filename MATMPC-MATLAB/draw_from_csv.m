%% load res in c

c_sim_res = readtable('data/sim_results.csv');

c_Ts = 0.025;
c_n_sim = length(c_sim_res.p);
c_t = 0:c_Ts:c_n_sim*c_Ts-c_Ts;

%%

figure(1);
subplot(321)
hold on
grid on
plot(c_t,c_sim_res.p,'--');
title('p');
subplot(322)
hold on
grid on
plot(c_t,c_sim_res.theta*180/pi,'--');
title('\theta');
subplot(323)
hold on
grid on
plot(c_t,c_sim_res.v,'--');
title('v');
subplot(324)
hold on
grid on
plot(c_t,c_sim_res.omega*180/pi,'--');
title('\omega');
subplot(3,2,[5 6]);
title('F');
hold on
grid on
stairs(c_t(2:end),c_sim_res.u(1:end-1),'--');
legend('matlab','c implementation')