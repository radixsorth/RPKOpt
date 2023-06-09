% initialization & data loading from files

global PERSIST_batchmode;
global PERSIST_fig_handle;


if(~isempty(PERSIST_batchmode) && PERSIST_batchmode==1)
    % batch mode ==> no graphical output
else
    % initialize the figure if not already initialized
    if(isempty(PERSIST_fig_handle))
        PERSIST_fig_handle = figure;
    else
        % only clear the figure and do not close the window
        clf(PERSIST_fig_handle);
    end
end

clearvars -except PERSIST*
clearvars -global -except PERSIST*

% global variables shared among the scripts
global plot_initial_guess;
global optimizer;
global optimization_instance epoch_within_instance beta_history;
global exp_count exp_names N_data rho_data final_time;
global beta Lambda lambda c0 N0;
global beta_size beta_sum_Serpent;
global seed_count seed_file shots rel_magnitude rand_magnitude;
global max_retries delta omega epsilon epsilon_a_ratio adaptivity grad_norm_mode options_primary adjoint options_adjoint learning_rate rel_learning_rate epoch epochs_to_go;
global loss_acc loss_acc_total;
global max_data_points_to_plot;

% -------------------------------------------------------------------------
% Model parameters

betas = 6;

switch betas
    
    case 6
        % --- 6-class version ---
        Lambda = 4.988e-5;
        lambda = [0.01334,0.03272,0.12081,0.30312,0.85096,2.85791];
        % initial guess for betas
        beta = [0.00027,0.00139,0.00134,0.00295,0.00122,0.00056];
        % intentionally wrong initial data
        %beta = [0.00237,0.00150,0.00089,0.00395,0.00162,0.00056];
        rand_magnitude = [0.00001, 0.0001, 0.0001, 0.0002, 0.0001, 0.00005];

    case 8
        % --- 8-class version ----
        Lambda = 4.202714e-5;
        % initial guess for betas
        lambda = [0.01247,0.02829,0.04252,0.13304,0.29247,0.66649,1.63478,3.55460];
        beta = [0.00031,0.0012,0.00076,0.00152,0.00263,0.00072,0.00066,0.00019];
        rand_magnitude = [0.00001, 0.0001, 0.0001, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001];

    otherwise
        disp('ERROR: Data for the given number of betas not available!')
        return
end

% the reactivity data are provided in the data files in beta_eff units
% (divided by "sum(beta)") which is determined by the Serpent code
% As sum(beta) is subject to optimization, we must use rho in absolute
% (dimensionless) form and thus pre-multiply it by the value below:
beta_sum_Serpent = 0.0077;

% -------------------------------------------------------------------------
% Common initial guess (seed) / randomization parameters

seed_count = 10;    % for reactor_seed.m
seed_file = 'beta_seed.txt';
shots = 50;        % for reactor_randomize.m
rel_magnitude = 0.2;

% ---------------------------
% number of delayed neutron classes ==> dimension of the optimization
% problem is beta_size, the number of equations is (beta_size+1)
beta_size = size(beta,2);
% ---------------------------

% -------------------------------------------------------------------------
% Automatic optimization parameters

plot_initial_guess = 1;         % if nonzero, the solutions for the initial guess will be added to the plots for comparison

optimizer = 'NesterovGD';        % valid choices: VanillaGD, ADAM, NesterovGD
epochs_to_go = 100;
adaptivity = 0;
adjoint = 1;                    % gradient computation: 0 = sensitivity analysis, 1 = adjoint method
grad_norm_mode = 1;             % 1 = gradient is normalized, 0 = gradient is divided by its initial norm upon calling reactor_optimize
learning_rate = 5e-6;
rel_learning_rate = ones(1,beta_size);
epsilon = ones(1,beta_size)*1e-6;   % fixed step for the finite difference approximation of gradient components
epsilon_a_ratio = 0.2;              % ... or the step is related to learning rate (see calc_gradient())
delta = 0.2;
omega = 0.8;
max_retries = 100;

% plot functions will resort to plotting every 2nd, 3rd etc. data point
% in the loss function & trajectories if the number of epochs is above this
% value
max_data_points_to_plot = 101;

% ODE solver options (different for fixed & adaptive learning rate)
% options = odeset('RelTol',1e-5,'AbsTol',1e-2,'Stats','on');
options_primary = { odeset('RelTol',1e-5,'AbsTol',1e-2) ...
            odeset('RelTol',1e-12,'AbsTol',1e-3) };
options_adjoint = { odeset('RelTol',1e-5,'AbsTol',1e-2, 'MaxStep', 0.1) ...
            odeset('RelTol',1e-12,'AbsTol',1e-3, 'MaxStep', 0.1) };

        
% -------------------------------------------------------------------------
% Experiments description & data files

% experiments to consider
experiments = [1 2 3];

% all available experiments listed here
exp_names = { 'Rod drop', 'Sine', 'Triangular' };
N_datafiles = { 'datafiles/pad_LCM_vyk_mov_avr.dat', 'datafiles/sin_LCM_vyk_mov_avr.dat', 'datafiles/troj100_LCM_vyk_mov_avr.dat' };
rho_datafiles = { 'datafiles/pad_LCM_rho.dat', 'datafiles/sin_LCM_rho.dat', 'datafiles/troj100_LCM_rho.dat' };


% ------------------------------------------------------------------------
% *** NO USER-ADJUSTABLE PARAMETERS BELOW THIS LINE ***
% ------------------------------------------------------------------------

% override the value of "experiments" (all other parameters can be modified
% AFTER the call to reactor_init, but this one needs to be known in advance
% as the data files are about to be loaded NOW!

global PERSIST_experiments
if(~isempty(PERSIST_experiments))
    experiments = PERSIST_experiments;
end

% only consider the subset of experiments given by the indices in the array
% 'experiments'
exp_names = exp_names(experiments);
N_datafiles = N_datafiles(experiments);
rho_datafiles = rho_datafiles(experiments);

exp_count = size(exp_names,2);

for j = 1:exp_count
    N_exp_data = readmatrix(N_datafiles{j});
    final_time{j} = N_exp_data(end,1);
    N_data{j} =  griddedInterpolant(N_exp_data(:,1), N_exp_data(:,2));
    rho_exp_data = readmatrix(rho_datafiles{j});
    % the reactivity data are provided in the data files in beta_eff units
    % (multiplied by beta_eff which is equal 
    rho_data{j} =  griddedInterpolant(rho_exp_data(:,1), rho_exp_data(:,2)*beta_sum_Serpent);
    % initial conditions
    N0{j} = N_data{j}(0);
    c0{j} = beta*N0{j}./(Lambda*lambda);
end

% initialization of the optimizer
epoch = 0;
optimization_instance = 1;
epoch_within_instance = 0;
beta_history{optimization_instance} = beta;

% valaue of the loss function (accumulated in order to be plotted)
% the initial value is set at 100%
for j = 1:exp_count
    loss_acc{j} = [100];
end
loss_acc_total = [100];
