global solved_once;

global epoch exp_count;
global epochs_to_go;
global optimizer;
global optimization_instance epoch epoch_within_instance beta_history;
global loss init_loss loss_total init_loss_total loss_acc loss_acc_total;

global grad_norm;

if(isempty(solved_once))
    reactor_solve;
end

% -------------------------------------------------------------------------

% MAIN LOOP

grad_norm = 0;

tic

for ep = 1:epochs_to_go
    fprintf('Processing epoch %d ... ',epoch+1);

    switch optimizer
        case 'VanillaGD'
            status = VanillaGD(ep);
        case 'ADAM'
            status = ADAM(ep);
        case 'NesterovGD'
            status = NesterovGD(ep);
        otherwise
            disp('Error: unknown optimizer');
            return
    end
    
    if status
        disp('Could not determine a suitable learning rate');
        break;
    end
    
    epoch = epoch + 1;
    epoch_within_instance = epoch_within_instance + 1;
    for j=1:exp_count
        loss_acc{j}(epoch+1) = 100 *  (loss(j) / init_loss(j));
    end
    loss_acc_total(epoch+1) = 100 *  (loss_total / init_loss_total);
    % log the evolution of beta
    beta_history{optimization_instance}(epoch_within_instance+1,:) = beta;

    fprintf('Partial relative losses: [ ');
    for j=1:exp_count
        fprintf('%g%% ', loss_acc{j}(end));
    end
    fprintf('], total loss %g%%.\n', loss_acc_total(end));
    
    reactor_plot;
    drawnow;

end

toc

disp('Done');



% -----------------

function gradient = calc_gradient_simple()
    % calculates the gradient of loss w.r.t. beta using sensitivity analysis
    % assumes that the the currently available solution already corresponds to the
    % curent value of beta
    global beta_size;
    global loss beta epsilon epsilon_a_ratio learning_rate rel_learning_rate;
    beta_orig = beta;
    l_orig = sum(loss);

    grad = zeros(1,beta_size);
    for i = 1:beta_size
        epsilon_a = rel_learning_rate(i) * learning_rate*epsilon_a_ratio;
        %epsilon_a = epsilon(i);
        beta(i) = beta_orig(i) + epsilon_a;
        reactor_solve;
        l=sum(loss);
        grad(i) = (l - l_orig)/epsilon_a;

        beta(i) = beta_orig(i);
    end
    
    gradient = grad;
end

% ---- ADJOINT BEGIN ----
function r = RHS_adjoint(j,T,t,pq)
    % NOTE: pq = [p q1 q2 q3 ... qX]
    global lambda beta_div_Lambda beta_sum Lambda rho_data;
    global beta_size;
    global exp_weights N_sim_interpolant N_data;
    r = zeros(beta_size+1,1);
    r(1) = ((rho_data{j}(T-t)-beta_sum)/Lambda)*pq(1) + beta_div_Lambda'*pq(2:beta_size+1) + exp_weights{j} * (N_data{j}(T-t) - N_sim_interpolant{j}(T-t));
    r(2:beta_size+1) = lambda'.*(pq(1)-pq(2:beta_size+1));
end
        
function gradient = calc_gradient_adjoint()
    % calculates the gradient by the adjoint method
    % assumes that the the currently available solution already corresponds to the
    % curent value of beta
    global beta_size;
    global exp_count final_time options_adjoint adaptivity;
    global time solution;
    global N_sim_interpolant Lambda lambda N0;

    % backward-in-time solution evaluation requires interpolation

    % solve the adjoint problems
    for j=1:exp_count
        N_sim_interpolant{j} = griddedInterpolant(time{j}, solution{j}(:,beta_size+1));
        [time_adj{j},pq{j}] = ode15s(@(t,pq) RHS_adjoint(j,final_time{j},t,pq),[0 final_time{j}],zeros(beta_size+1,1), options_adjoint{adaptivity+1});
    end

    % We take advantage of the fact that the solution interpolant is already available
    % and we perform the integration in the reverse time direction (along
    % the solution of the transformed adjoint problem)
    % All components of the loss function gradient w.r.t. beta are calculated at once.
    grad = zeros(1,beta_size);
    for j=1:exp_count
        grad = grad + ( trapz(time_adj{j},(pq{j}(:,1)-pq{j}(:,2:beta_size+1)).*N_sim_interpolant{j}(final_time{j}-time_adj{j}),1) ) / Lambda;
        % contribution of the initial condition
        grad = grad - N0{j}./(Lambda*lambda).*pq{j}(end,2:beta_size+1);
    end
    
    gradient = grad;
end
% ---- ADJOINT END ----
        
function gradient = calc_gradient()
    global adjoint grad_norm grad_norm_mode;
    if(adjoint)
        gradient = calc_gradient_adjoint();
    else
        gradient = calc_gradient_simple();
    end
    
    % normalize the gradient, using its current or initial (upon calling this script) value
    % (depends on the value of 'grad_norm_mode')
    if(grad_norm_mode || grad_norm==0)
        grad_norm=norm(gradient);
    end
    gradient = gradient / grad_norm;
end

% -------------------------------------------------------------------------
% Vanilla GD
% -------------------------------------------------------------------------

function stat = VanillaGD(ep) % one (successful) step of GD
    
    global adaptivity max_retries learning_rate rel_learning_rate delta omega;
    global beta time solution;
    persistent gradient;

    if(~adaptivity || ep==1)
            gradient = calc_gradient();
    end

    if(adaptivity)
        % store the current state (not used if adaptivity is off)
        prev_beta = beta;
        prev_gradient = gradient;
        prev_time = time;
        prev_solution = solution;
    end

    for retries = 0:max_retries

        % proceed a step further with the current learning rate
        beta = beta - learning_rate * rel_learning_rate.*gradient;
        reactor_solve;

        if(adaptivity)
            fprintf('Retries: %d, learning rate = %g\n', retries, learning_rate);
            gradient = calc_gradient();
            e = 0.5 * (1-dot(gradient, prev_gradient)/(norm(gradient)*norm(prev_gradient)));

            % adjust the learning rate
            if(e>0.001)
                learning_rate = (delta/e)^0.2 * omega * learning_rate;
            else
                % if e is too close to zero
                learning_rate = 2*learning_rate;
            end
        end

        if(~adaptivity || e<delta)
            % successful step ==> move on to the next epoch
            stat = 0;
            return;
        end

        % go back and try again
        beta = prev_beta;
        gradient = prev_gradient;
        time = prev_time;
        solution = prev_solution;
    end

    stat = 1;
end

% -------------------------------------------------------------------------
% ADAM
% -------------------------------------------------------------------------

function stat = ADAM(ep) % one (successful) step of ADAM GD
    
    global learning_rate;
    global beta;
    persistent A F;
    
    % ADAM parameters
    rho = 0.999;
    rho_f = 0.9;
    eps_reg = 1e-8;

    if(ep==1)
            % initialize ADAM
            A = zeros(size(beta));
            F = A;
    end

    gradient = calc_gradient();
    A = rho*A + (1-rho)*gradient.^2;
    F = rho_f*F + (1-rho_f)*gradient;
    alpha = learning_rate * sqrt(1-rho^ep) / (1-rho_f^ep);
    
    beta = beta - alpha./sqrt(A+eps_reg).*F;
    
    reactor_solve;

    % successful step ==> move on to the next epoch
    stat = 0;
end

% -------------------------------------------------------------------------
% Nesterov GD (optionally with gradient evaluation averaging) and with adaptivity tailored
% to the Nesterov algorithm
% -------------------------------------------------------------------------

function [gradient, e] = calc_gradient_nesterov(ep)

    global beta time solution learning_rate rel_learning_rate grad_norm grad_norm_mode;
    persistent gradient1;

    % as always, the solution with the current beta is expected to be
    % readily available
    curr_beta = beta;
    curr_time = time;
    curr_solution = solution;
    
    % Nesterov momentum parameter (relative to the learning rate!)
    momentum = 1.0;
    % additional 'look-ahead distance' parameter (relative to the learning rate!)
    % (a parameter for which a value 0.5 reminds of the 2nd order midpoint Runge-Kutta method)
    % In the original Nesterov momentum GD algorithm, it is equal to
    % momentum.
    look_ahead = 0.5;
    
    % force the modified Nesterov momentum GD where the gradient from the
    % previous iteration is replaced by the gradient at the current
    % beta (requires 2 gradient evaluations per epoch)
    % For this algorithm, momentum=1.0 is the best choice as it most
    % efficiently cancels oscillations (which is also related to using
    % normalized gradients all the time - see 'grad_norm_mode' in reactor_init)
    ep = 1;

    if(ep==1)
        % initialize Nesterov momentum GD
        gradient1 = calc_gradient();
    end
    
    % look ahead
    beta = beta - look_ahead*learning_rate * rel_learning_rate.*gradient1;
    reactor_solve;
    gradient2 =  calc_gradient();

    % calculate the deviation of the two successive gradients to allow
    % learning rate adaptivity
    e = 0.5 * (1-dot(gradient1, gradient2));

    % Nesterov GD with averaging (always normalizing the gradient)
    gradient = momentum * gradient1 + gradient2;
    if(grad_norm_mode || grad_norm==0)
        gradient = gradient / norm(gradient);
    end
    
    % next time, we will re-use the current direction
    gradient1 = gradient;
    
    % revert to the state before gradient evaluation
    beta = curr_beta;
    time = curr_time;
    solution = curr_solution;
end

function stat = NesterovGD(ep) % one (successful) step of GD
    
    global adaptivity max_retries learning_rate rel_learning_rate delta omega;
    global beta;

    for retries = 0:max_retries

        [gradient, e] = calc_gradient_nesterov(ep);

        if(adaptivity)
            fprintf('Retries: %d, learning rate = %g\n', retries, learning_rate);
            % adjust the learning rate
            if(e>0.001)
                learning_rate = (delta/e)^0.2 * omega * learning_rate;
            else
                % if e is too close to zero
                learning_rate = 2*learning_rate;
            end

            if(e>delta)
                % deviation too large; repeat with a smaller learning rate
                continue;
            end
        end
         
        % proceed a step further with the current learning rate
        beta = beta - learning_rate * rel_learning_rate.*gradient;

        % this is actually not needed in the classical Nesterov GD, but we do that
        % 1) to evaluate & plot the loss at the current position
        % 2) to always finish with the solution evaluated at current beta (for compatibility with the other methods)
        % 3) to allow the modified Nesterov GD algorithm (see above)
        reactor_solve;

        % successful step ==> move on to the next epoch
        stat = 0;
        return;
    end

    stat = 1;
end
