% solution of the system with the current settings of beta

global beta Lambda lambda N0 time solution;
global beta_size;
global beta_sum beta_div_Lambda;
global exp_count N_data final_time;
global loss init_loss loss_total init_loss_total adaptivity options_primary;
global init_time init_solution;
global solved_once;

beta_sum = sum(beta);
beta_div_Lambda = beta' / Lambda;

% perform the computation 
for j=1:exp_count
    % calculate the initial condition that depends on the current beta
    c0{j} = beta*N0{j}./(Lambda*lambda);
    % solve the ODE system
    [time{j},solution{j}] = ode15s(@(t,y) RHS_primary(j,t,y),[0 final_time{j}],[c0{j} N0{j}], options_primary{adaptivity+1});
end

% update the losses
for j=1:exp_count
        loss(j) = 0.5*trapz(time{j}, (solution{j}(:,beta_size+1)-N_data{j}(time{j})).^2);
end
loss_total = sum(loss);

% initial solution - store the 
if(isempty(solved_once))
    init_loss = loss;
    init_loss_total = loss_total;
    
    % regularization for calculation of the relative value of loss
    init_loss(init_loss<1e-6) = 1e-6;
    if(init_loss_total<1e-6)
        init_loss_total = 1e-6;
    end
    
    % save the solution for the initial guess in order to be able to include
    % it in the plots (if plot_initial_guess==1)
    init_time = time;
    init_solution = solution;
    
    solved_once = 1;
end


% -------------------------

function r = RHS_primary(j,t,y)
            global lambda Lambda beta_sum beta_div_Lambda rho_data;
            global beta_size;
            % j ... experiment number
            % y(t) = [c1(t), c2(t), c3(t), .... , cX(t), N(t)]
            r = zeros(beta_size+1,1);
            r(1:beta_size) = beta_div_Lambda*y(beta_size+1) - lambda'.*y(1:beta_size);
            r(beta_size+1) = ((rho_data{j}(t)-beta_sum)/Lambda)*y(beta_size+1) + dot(lambda,y(1:beta_size));
end
