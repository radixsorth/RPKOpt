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
global exp_count exp_names exp_weights N_data rho_data final_time;
global beta beta_uncertainties drift_penalization Lambda lambda N0;
global beta_init beta_size beta_sum_Serpent;
global seed_count seed_file shots rel_magnitude rand_magnitude;
global max_retries delta omega epsilon epsilon_a_ratio adaptivity grad_norm_mode options_primary adjoint options_adjoint learning_rate rel_learning_rate epoch epochs_to_go;
global loss_acc loss_acc_total;
global max_data_points_to_plot;

% -------------------------------------------------------------------------
% Model parameters

opt_case = 'C15_6';
raw_data = 0;           % 0 = moving-averaged power output data, 1 = raw data

switch opt_case
    
    case 'C15_6'
        % --- 6-class version ---
        core_configuration = 'C15';
        Lambda = 4.19e-5;
        lambda = [0.01334,0.03272,0.12081,0.30312,0.85096,2.85791];
        % initial guess for betas
        beta = 1e-4 * [2.72 13.9 13.4 29.8 12.4 5.18];
        % intentionally wrong initial data
        %beta = [0.00237,0.00150,0.00089,0.00395,0.00162,0.00056];
        beta_uncertainties = 1e-6 * [1.75 3.96 3.98 5.83 3.84 2.45];
        rand_magnitude = [0.00001, 0.0001, 0.0001, 0.0002, 0.0001, 0.00005];

    case 'C15_8'
        % TODO
        % --- 8-class version ----
        core_configuration = 'C15';
        Lambda = 4.18e-5;
        % initial guess for betas
        lambda = [0.01247,0.02829,0.04252,0.13304,0.29247,0.66649,1.63478,3.55460];
        % orig. data for 8 classes
        %beta = [0.00031,0.0012,0.00076,0.00152,0.00263,0.00072,0.00066,0.00019];
        % alternate data for 8 classes
        beta = 1e-4 * [2.71 11.8 7.5 15.7 25.9 7.4 6.54 1.91];
        beta_uncertainties = 1e-6 * [1.74 3.68 2.89 4.32 5.61 2.90 2.70 1.47];
        rand_magnitude = [0.00001, 0.0001, 0.0001, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001];

    case 'C12-C_6'
        % --- 6-class version ---
        core_configuration = 'C12-C';
        Lambda = 4.93e-5;
        lambda = [0.01334,0.03272,0.12081,0.30312,0.85096,2.85791];
        % initial guess for betas
        beta = 1e-4 * [2.69 13.8 13.3 29.4 12.2 5.14];
        % intentionally wrong initial data
        %beta = [0.00237,0.00150,0.00089,0.00395,0.00162,0.00056];
        beta_uncertainties = 1e-6 * [2.80 6.21 6.43 9.11 6.15 3.87];
        rand_magnitude = [0.00001, 0.0001, 0.0001, 0.0002, 0.0001, 0.00005];

    case 'C12-C_8'
        % TODO
        % --- 8-class version ----
        core_configuration = 'C12-C';
        Lambda = 4.93e-5;
        % initial guess for betas
        lambda = [0.01247,0.02829,0.04252,0.13304,0.29247,0.66649,1.63478,3.55460];
        % orig. data for 8 classes
        %beta = [0.00031,0.0012,0.00076,0.00152,0.00263,0.00072,0.00066,0.00019];
        % alternate data for 8 classes
        beta = 1e-4 * [2.66 11.7 7.45 15.5 25.5 7.19 6.45 1.87];
        beta_uncertainties = 1e-6 * [2.76 5.82 4.44 6.49 8.7 4.65 4.32 2.29];
        rand_magnitude = [0.00001, 0.0001, 0.0001, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001];

    otherwise
        disp('ERROR: Data for the given case not available!')
        return
end

% -------------------------------------------------------------------------
% Experiments description & data files

% weight on the regularization term (gamma)
drift_penalization = 1e4;

% the reactivity data are provided in the data files in beta_eff units
% (divided by "sum(beta)") which is determined by the Serpent code
% As sum(beta) is subject to optimization, we must use rho in absolute
% (dimensionless) form and thus pre-multiply it by the value below:
beta_sum_Serpent = 0.0077;


% experiments to consider
experiments = [1 2 3];

% all available experiments listed here
exp_names = { 'Rod drop', 'Sine', 'Triangular' };
datafiles = {'drop.dat', 'sin.dat', 'triang.dat'};
% weights of the individual experiments in the loss functional
% (below, they are further adjusted to account for loss per unit time of each experiment)
exp_weights = {1, 1, 1};

raw_prefix = {'', 'raw_'};
% path to power output datafiles
N_path = ['datafiles/' core_configuration '/' raw_prefix{raw_data+1} 'power'];
% path to reactivity datafiles
rho_path = ['datafiles/' core_configuration '/reactivity'];

% scaling factor between raw measurements (milliamperes) and the correct
% power output units (impulses per second)
mA_to_impps = 1e5/105.0;

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
learning_rate = 2.5e-5;
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
options_primary = { odeset('RelTol',1e-5,'AbsTol',1e+1) ...
            odeset('RelTol',1e-9,'AbsTol',1e-3) };
options_adjoint = { odeset('RelTol',1e-5,'AbsTol',1e-3, 'MaxStep', 0.1) ...
            odeset('RelTol',1e-9,'AbsTol',1e-7, 'MaxStep', 0.1) };


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
datafiles = datafiles(experiments);
exp_weights = exp_weights(experiments);

exp_count = size(exp_names,2);

for j = 1:exp_count
    N_exp_data = readmatrix([N_path '/' datafiles{j}]);
    final_time{j} = N_exp_data(end,1);
    N_data{j} =  griddedInterpolant(N_exp_data(:,1), N_exp_data(:,2)*mA_to_impps);
    rho_exp_data = readmatrix([rho_path '/' datafiles{j}]);
    % the reactivity data are provided in the data files in beta_eff units
    % (multiplied by beta_eff). We convert it to dimensionless units.
    rho_data{j} =  griddedInterpolant(rho_exp_data(:,1), rho_exp_data(:,2)*beta_sum_Serpent);
    % initial conditions
    N0{j} = N_data{j}(0);
    % adjust experiment weights by the inverse proportionality to final time
    % (relate loss to one second of time for each experiment)
    exp_weights{j} = exp_weights{j} / final_time{j};
end

% initialization of the optimizer
epoch = 0;
optimization_instance = 1;
epoch_within_instance = 0;
beta_history{optimization_instance} = beta;
beta_init = beta;

% valaue of the loss function (accumulated in order to be plotted)
% the initial value is set at 100%
for j = 1:exp_count
    loss_acc{j} = [100];
end
loss_acc_total = [100];
