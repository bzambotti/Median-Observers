function plot_results(x,x_hat,e_norm,...
                      y_clean,y_corrupted,...
                      disturbance,...
                      selected_delay,...
                      innovation_raw,...
                      innovation_filtered)


n = size(x,1);
p = size(y_clean,1);

sim_steps = size(y_clean,2);

time = 0:sim_steps;
tk   = 1:sim_steps;

figure('Color','w',...
       'Position',[100 50 1300 950]);

subplot(5,1,1)

semilogy(time,e_norm,...
    'LineWidth',2)

grid on
box on

xlabel('k')
ylabel('||e_k||')

title('Estimation error')

subplot(5,1,2)

hold on

colors=lines(n);

for i=1:n

    plot(time,...
        x(i,:),...
        '-',...
        'Color',colors(i,:),...
        'LineWidth',2);

    plot(time,...
        x_hat(i,:),...
        '--',...
        'Color',colors(i,:),...
        'LineWidth',2);

end

grid on
box on

xlabel('k')

title('True and estimated states')

L={};

for i=1:n
    L{end+1}=sprintf('x_%d',i);
    L{end+1}=sprintf('\\hat{x}_%d',i);
end

legend(L,'Location','eastoutside')


subplot(5,1,3)

hold on

for j=1:p

    plot(tk,...
        y_clean(j,:),...
        'k',...
        'LineWidth',2);

    plot(tk,...
        y_corrupted(j,:),...
        'rx',...
        'MarkerSize',7,...
        'LineWidth',1.5);

end

grid on
box on

xlabel('k')

title('Outputs')

legend('True output','Corrupted output')


subplot(5,1,4)

hold on

if p==1

    stairs(tk,...
           selected_delay,...
           'LineWidth',2);

else

    colors=lines(p);

    for j=1:p

        stairs(tk,...
               selected_delay(j,:),...
               'Color',colors(j,:),...
               'LineWidth',2);

    end

end

ylim([-0.5 max(selected_delay(:))+0.5])

yticks(0:max(selected_delay(:)))

grid on
box on

xlabel('k')

ylabel('Delay')

title('Delay selected by the median observer')


subplot(5,1,5)

hold on

plot(tk,...
     innovation_raw(1,:),...
     '--',...
     'LineWidth',1.2);

plot(tk,...
     innovation_filtered(1,:),...
     'LineWidth',2);

if nargin>=6

    yyaxis right

    stem(tk,...
         disturbance(1,:),...
         'filled');

    ylabel('Outlier')

end

grid on
box on

xlabel('k')

title('Innovation before and after median filtering')

legend('Raw innovation',...
       'Median innovation',...
       'Outlier')

sgtitle('Median Observer Simulation')