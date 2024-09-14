% A study of the influence of ODE solver accuracy and gradient computation
% method

global PERSIST_batchmode
% suppress graphics output
PERSIST_batchmode=1;

% variables set in reactor_init.m and overriden in this study need to be
% declared global
global optimizer adjoint options_primary options_adjoint learning_rate PERSIST_experiments

PERSIST_basedir = 'Direct-vs-adjoint/';

% solver
PERSIST_optimizer = 'VanillaGD';

% accuracy
PERSIST_RelTol = [ 1e-5, 1e-7, 1e-9 ];  % both for primary & adjoint
PERSIST_AbsTol = [ 1e+1, 1e-1, 1e-3 ];  % for primary only; value for adjoint is 10^4 times smaller
PERSIST_acc_name = { 'low', 'medium', 'high' };

% gradient computation method
PERSIST_adjoint = [ 0, 1 ];
PERSIST_grad_method_name = { 'direct', 'adjoint' };

% experiments
PERSIST_experiment_set = { [1], [1 2 3] };
PERSIST_exp_name = { 'drop', 'all' };
PERSIST_learning_rate = [ 5e-6, 1e-6 ];

for PERSIST_acc =1:3
    for PERSIST_method = 1:2
        for PERSIST_exp = 1:2

            % initialize case
            clc
            PERSIST_experiments = PERSIST_experiment_set{PERSIST_exp};  % override "experiments" in reactor_init
            % this clears all variables not named PERSIST_*
            reactor_init
            case_name = [ PERSIST_exp_name{PERSIST_exp} '_' PERSIST_optimizer '_' PERSIST_grad_method_name{PERSIST_method} '_' PERSIST_acc_name{PERSIST_acc} ];

            % override all remaining case parameters
            optimizer = PERSIST_optimizer;
            adjoint = PERSIST_adjoint(PERSIST_method);
            options_primary = { odeset('RelTol', PERSIST_RelTol(PERSIST_acc) ,'AbsTol', PERSIST_AbsTol(PERSIST_acc)) ...
            odeset('RelTol', PERSIST_RelTol(PERSIST_acc)/10 ,'AbsTol', PERSIST_AbsTol(PERSIST_acc)/10) };
            options_adjoint = { odeset('RelTol', PERSIST_RelTol(PERSIST_acc) ,'AbsTol', PERSIST_AbsTol(PERSIST_acc)/1e4, 'MaxStep', 0.1) ...
            odeset('RelTol', PERSIST_RelTol(PERSIST_acc)/10 ,'AbsTol', PERSIST_AbsTol(PERSIST_acc)/10/1e4, 'MaxStep', 0.1) };
            learning_rate = PERSIST_learning_rate(PERSIST_exp);

            % turn on logging to file
            diary([ PERSIST_basedir case_name '.log']);

            % run the optimization case
            disp([ 'Case:' case_name ]);
            reactor_optimize

            % turn logging off
            diary off

            % save the workspace
            save([ PERSIST_basedir case_name '.mat']);

        end
    end
end

disp('All cases finished.');
