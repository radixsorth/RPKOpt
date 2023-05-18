% run optimizations from randomized starting points
% may be re-run multiple times

% If run immediately after reactor_initialize, the seed_file is examined
% and if it exists, starting points are read from that file.
% Otherwise, random seeds are generated and stored in seed_file.

% Each subsequent call to reactor_seed adds more random seeds that are not
% stored to file.

paths=figure;

global epoch optimization_instance epoch_within_instance beta_history;
global beta beta_size beta_init;
global seed_count seed_file rel_magnitude rand_magnitude;

reseed = 1;
if(epoch==0)
    beta_init = beta;

    if(isfile(seed_file))
        betas = readmatrix(seed_file);
        seed_count = size(betas,1);
        beta_sz = size(betas,2);
        if(beta_sz == beta_size)
            disp('Initial seeds for beta have beeen loaded from file.');
            reseed = 0;
        else
            disp('Warning: beta seed data does not match the dimension of the problem! Will re-seed.');
        end
    end
end

if(reseed)
    fprintf('Generating %d random seeds.\n', seed_count);
    % generate random seeds around the initial beta
    betas = zeros(seed_count, beta_size);
    for i=1:seed_count
        % the very first pass is with the original beta
        betas(i,:) = beta_init + (epoch>0 || i>1)*rel_magnitude*rand_magnitude.*(2*rand(1,beta_size)-1);
    end
    
    if(epoch==0)
        % write the random initial seeds to file for further re-use
        fprintf('Saving random initial points to: %s\n', seed_file);
        writematrix(betas,seed_file);
    end
end

for i=1:seed_count

    beta = betas(i,:);
    % prepare the solution with the starting beta
    reactor_solve;

    if(epoch>0)
        % the initial value of beta in each optimization instance
        % is treated as a new epoch. The only exception occurs when
        % reactor_seed is run immediately after reactor_initialize
 
        % reset epoch counter within the current optimization instance
        % (only used to store beta history)
        epoch_within_instance = 0;
        
        epoch = epoch + 1;
        for j=1:exp_count
            loss_acc{j}(epoch+1) = 100 *  (loss(j) / init_loss(j));
        end
        loss_acc_total(epoch+1) = 100 *  (loss_total / init_loss_total);

        % log the evolution of beta
        optimization_instance = optimization_instance + 1;
        beta_history{optimization_instance} = beta;
        fprintf('---- Optimization instance %d ----\n', optimization_instance);
    end

    fprintf('--> Seed %d (epoch %d): Partial relative losses: [ ', optimization_instance, epoch);
    for j=1:exp_count
        fprintf('%g%% ', loss_acc{j}(end));
    end
    fprintf('], total loss %g%%.\n', loss_acc_total(end));

    reactor_plot;
    drawnow;

    % run the GD optimizer
    reactor_optimize;

    fprintf('Completed, plotting optimization path.\n');
    figure(paths);
    reactor_plot_paths;
    drawnow;
end
