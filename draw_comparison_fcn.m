function [fig1, fig2] = draw_comparison_fcn(time,state_sim,one_step_pred,controls_MPC,color)
%DRAW_COMPARISON_FCN Summary of this function goes here
%   Detailed explanation goes here

fig1 = figure(1);
subplot(311)
hold on
grid on
plot(time,state_sim(:,1),color);
ylabel('p [m]');
subplot(312)
hold on
grid on
plot(time,state_sim(:,2)*180/pi,color);
ylabel('\theta [deg]')
% subplot(313)
% hold on
% grid on
% plot(time,state_sim(:,3));
% legend('sim','one-step pred')
% title('v');
% subplot(324)
% hold on
% grid on
% plot(time,state_sim(:,4)*180/pi);
% plot([time(1:end-1)', time(2:end)']',[state_sim(1:end-1,4),one_step_pred(:,4)]'*180/pi,'--')
% legend('sim','one-step pred')
% title('\omega');
subplot(3,1,3);
ylabel('F [N]');
hold on
grid on
stairs(time(1:end-1),controls_MPC(:,1),color);
xlabel('Time [s]')


fig2 = figure(2);
% subplot(221)
% hold on
% grid on
% plot(time(2:end),state_sim(2:end,1)-one_step_pred(:,1),color);
% % title('difference sim vs pred - p');
% subplot(222)
% hold on
% grid on
% plot(time(2:end),state_sim(2:end,2)-one_step_pred(:,2),color);
% title('difference sim vs pred - \theta');
subplot(211)
hold on
grid on
plot(time(2:end),abs(state_sim(2:end,3)-one_step_pred(:,3)),color);
ylabel('|\Delta v| [m/s]')
% title('difference sim vs pred - v');
subplot(212)
hold on
grid on
plot(time(2:end),abs(state_sim(2:end,4)-one_step_pred(:,4)),color);
ylabel('|\Delta \omega| [rad/s]')
xlabel('Time [s]')
% title('difference sim vs pred - \omega');

done = true;
end

