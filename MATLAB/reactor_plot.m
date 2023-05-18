global PERSIST_batchmode;

% batch mode ==> no graphical output
if(~isempty(PERSIST_batchmode))
    if(PERSIST_batchmode==1)
        return;
    end
end

exp_plot_resolution = 200;

global exp_count exp_names N_data;
global beta_size;
global time solution;
global epoch loss_acc loss_acc_total;
global ax tab ax_loss tab_loss;
global max_data_points_to_plot;
global solved_once plot_initial_guess;
global PERSIST_fig_handle;

if(isempty(solved_once))
    reactor_solve;
end

figure_exists = 0;
if(ishandle(ax_loss))
    if(isvalid(ax_loss))
        % figure already initialized and it still exists
        figure_exists = 1;
    end
end

% (re)create the figure with as many tabs as there are experiments
if(~figure_exists)
    % activate the figure dedicated for us
    figure(PERSIST_fig_handle);
    for j=1:exp_count
        tab{j} = uitab('Title',exp_names{j});
        ax{j} = axes(tab{j});
    end
    tab_loss = uitab('Title','Loss function');
    ax_loss = axes(tab_loss);
end


for j=1:exp_count
    hold(ax{j},'off');
    aux_time = linspace(0,final_time{j},exp_plot_resolution);
   
    plot(ax{j}, aux_time, N_data{j}(aux_time), 'LineWidth', 2);
    hold(ax{j},'on');
    plot(ax{j}, time{j}, solution{j}(:,beta_size+1), 'LineWidth', 2);
    if(plot_initial_guess)
        plot(ax{j}, init_time{j}, init_solution{j}(:,beta_size+1), '--', 'Color', '#888888', 'LineWidth', 1);
    end

    title(ax{j},[exp_names{j} ': Power output response'])
    xlabel(ax{j}, 'time (s)')
    ylabel(ax{j}, 'power output (imp/s)')
    if(plot_initial_guess)
        legend(ax{j},'experiment','simulation','simulation (initial guess)')
    else
        legend(ax{j},'experiment','simulation')
    end
end

% skip some data points in order to plot at most 'max_data_points_to_plot'
every_nth = ceil((epoch+1) / max_data_points_to_plot);

hold(ax_loss, 'off');
plot(ax_loss, 0:every_nth:epoch, loss_acc_total(1:every_nth:end), '-o', 'LineWidth', 2);
loss_legend_entries{1} = 'Total loss:';
hold(ax_loss, 'on');
for j=1:exp_count
       plot(ax_loss, 0:every_nth:epoch, loss_acc{j}(1:every_nth:end), '--o', 'LineWidth', 2);
       loss_legend_entries{j+1} = [ exp_names{j} ' loss' ];
end
xlabel(ax_loss, 'Epoch');
ylabel(ax_loss, 'Relative loss w.r.t. the initial beta [%]');

legend(ax_loss, loss_legend_entries);

        
