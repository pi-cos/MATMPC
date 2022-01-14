close all

%%

load sim_white_box_nom

figure(2);
subplot(221)
hold on
grid on
plot(time(2:end),state_sim(2:end,1)-one_step_pred(:,1));
title('difference sim vs pred - p');
subplot(222)
hold on
grid on
plot(time(2:end),state_sim(2:end,2)-one_step_pred(:,2));
title('difference sim vs pred - \theta');
subplot(223)
hold on
grid on
plot(time(2:end),state_sim(2:end,3)-one_step_pred(:,3));
title('difference sim vs pred - v');
subplot(224)
hold on
grid on
plot(time(2:end),state_sim(2:end,4)-one_step_pred(:,4));
title('difference sim vs pred - \omega');

%%

load sim_cont_grey_box

figure(2);
subplot(221)
hold on
grid on
plot(time(2:end),state_sim(2:end,1)-one_step_pred(:,1));
title('difference sim vs pred - p');
subplot(222)
hold on
grid on
plot(time(2:end),state_sim(2:end,2)-one_step_pred(:,2));
title('difference sim vs pred - \theta');
subplot(223)
hold on
grid on
plot(time(2:end),state_sim(2:end,3)-one_step_pred(:,3));
title('difference sim vs pred - v');
subplot(224)
hold on
grid on
plot(time(2:end),state_sim(2:end,4)-one_step_pred(:,4));
title('difference sim vs pred - \omega');

%%

load sim_cont_black_box

figure(2);
subplot(221)
hold on
grid on
plot(time(2:end),state_sim(2:end,1)-one_step_pred(:,1));
title('difference sim vs pred - p');
subplot(222)
hold on
grid on
plot(time(2:end),state_sim(2:end,2)-one_step_pred(:,2));
title('difference sim vs pred - \theta');
subplot(223)
hold on
grid on
plot(time(2:end),state_sim(2:end,3)-one_step_pred(:,3));
title('difference sim vs pred - v');
subplot(224)
hold on
grid on
plot(time(2:end),state_sim(2:end,4)-one_step_pred(:,4));
title('difference sim vs pred - \omega');

%%

legend('white box','grey box', 'black_box')