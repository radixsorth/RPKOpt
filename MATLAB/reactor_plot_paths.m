global PERSIST_batchmode;

% batch mode ==> no graphical output
if(~isempty(PERSIST_batchmode))
    if(PERSIST_batchmode==1)
        return;
    end
end

global PERSIST_fig_handle;
global beta beta_init beta_uncertainties beta_history epochs_to_go;

% avoid overwriting the content of the main figure
if(gcf==PERSIST_fig_handle)
    figure;
end

clf;

% plot the trajectories of the individual beta components
% we assume that the size of beta is even
subplot_count = size(beta,2);
 for pl=1:subplot_count
    subplot(subplot_count/2,2,pl);
    hold on;
    
    records = size(beta_history,2);

    % shade out the uncertainty range
    rectangle('Position',[0 beta_init(pl)-beta_uncertainties(pl) epochs_to_go 2*beta_uncertainties(pl)],'LineStyle','--','FaceColor', [0 0 0 0.07]);

    % plot the paths for all available optimization instances
    for r=1:records
        final_epoch = size(beta_history{r},1)-1;

        % skip some data points in order to plot at most 'max_data_points_to_plot'
        %every_nth = ceil((final_epoch+1) / max_data_points_to_plot);
        %plot(0:every_nth:final_epoch, beta_history{r}(1:every_nth:end,pl),'-o');

        % plot curves only
        plot(0:final_epoch, beta_history{r}(:,pl),'-','LineWidth',2);
        % disable auto x-range adjustment upon window resize
        xlim([0 final_epoch]);

    end
    xlabel('Epoch');
    ylabel(sprintf('$\\beta_{\\mathrm{eff},%d}$',pl),'Interpreter','latex');
 end
