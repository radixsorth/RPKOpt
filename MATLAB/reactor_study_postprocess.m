% Postprocess a study of the influence of ODE solver accuracy and gradient computation
% method

clear;

PERSIST_basedir = 'Direct-vs-adjoint/';

% solver
PERSIST_optimizer = 'VanillaGD';

% accuracy
PERSIST_acc_name = { 'low', 'medium', 'high' };
PERSIST_acc_color = { '#0072BD', '#D95319', '#77AC30' };

% gradient computation method
PERSIST_adjoint = [ 0, 1 ];
PERSIST_grad_method_name = { 'direct', 'adjoint' };
PERSIST_grad_method_linestyle = { '-', '--' };

% experiments
PERSIST_exp_name = { 'drop', 'all' };
PERSIST_learning_rate = [ 5e-6, 1e-6 ];

i=0;

% load data
for PERSIST_acc = 1:3
    for PERSIST_method = 1:2
        % it's best to evaluate per experiment (either "1" or "2")
        for PERSIST_exp = 2
            
            i=i+1;

            case_name = [ PERSIST_exp_name{PERSIST_exp} '_' PERSIST_optimizer '_' PERSIST_grad_method_name{PERSIST_method} '_' PERSIST_acc_name{PERSIST_acc} ];
            result{i}=load([ PERSIST_basedir case_name '.mat'],'beta_history','loss_acc_total');
            % customize the plot properties for each combination of
            % PERSIST_acc and PERSIST_method
            legend_entries{i}=[ PERSIST_grad_method_name{PERSIST_method} '_' PERSIST_acc_name{PERSIST_acc} ];
            linestyle{i} = PERSIST_grad_method_linestyle{PERSIST_method};
            color{i} = PERSIST_acc_color{PERSIST_acc};

        end
    end
end

study_count = i;

% -------------------------------------------------------------------------

% plot total loss

max_data_points_to_plot = 101;
epochs=size(result{1}.loss_acc_total,2);
figure

% skip some data points in order to plot at most 'max_data_points_to_plot'
every_nth = ceil((epochs+1) / max_data_points_to_plot);

hold on;
for i=1:study_count
    plot(0:every_nth:epochs, result{i}.loss_acc_total(1:every_nth:end), linestyle{i}, 'Color', color{i}, 'LineWidth', 2);
end
xlabel('Epoch');
ylabel('Relative loss w.r.t. the initial guess [%]','Interpreter','none');
legend(legend_entries,'Interpreter','none');

% plot the trajectories of the individual beta components
% we assume that the size of beta is even
figure
subplot_count = size(result{1}.beta_history{1},2);
 for pl=1:subplot_count
    subplot(subplot_count/2,2,pl);
    hold on;
    
    for i=1:study_count
        final_epoch = size(result{i}.beta_history{1},1)-1;

        % plot curves only
        plot(0:final_epoch, result{i}.beta_history{1}(:,pl),linestyle{i}, 'Color', color{i});
    end
    xlabel('Epoch');
    ylabel(sprintf('$\\beta_{\\mathrm{eff},%d}$',pl),'Interpreter','latex');
 end
legend(legend_entries,'Interpreter','none');
