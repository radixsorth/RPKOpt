% plots the evolution of the total loss function as multiple curves for
% each optimization instance

global PERSIST_fig_handle;
global max_data_points_to_plot ;
global optimization_instance epochs_to_go epoch loss_acc_total;

% to proceed, all optimization instances must have finished epochs_to_go
if(optimization_instance<2  ||  (epoch+1)/optimization_instance~=(epochs_to_go+1))
    disp('Not ready to plot loss progress over multiple optimization instances.');
    return;
end

% avoid overwriting the content of the main figure
if(gcf==PERSIST_fig_handle)
    figure;
end

clf;

% skip some data points in order to plot at most 'max_data_points_to_plot'
every_nth = ceil((epochs_to_go+1) / max_data_points_to_plot);

hold on;
for i=1:optimization_instance
    plot(0:every_nth:epochs_to_go, loss_acc_total((i-1)*(epochs_to_go+1)+1 : every_nth : i*(epochs_to_go+1)), '-o', 'LineWidth', 2);
end
xlabel('Epoch within instance');
ylabel('Relative loss w.r.t. the initial beta [%]');
