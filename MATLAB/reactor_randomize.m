% optimize by random shooting and moving in the direction of any
% improvement

global solved_once;
global shots rel_magnitude rand_magnitude;
global optimization_instance epoch_within_instance beta_history;
global exp_count beta loss loss_total init_loss_total loss_acc loss_acc_total epoch;

if(isempty(solved_once))
    reactor_solve;
end

best_beta = beta;
best_loss_total = loss_total;

for i=1:shots
    fprintf('Shoot %d ... ', i);
    beta = best_beta + rel_magnitude*rand_magnitude.*(2*rand(size(beta))-1);
    reactor_solve;
    fprintf('current loss = %g%%, best loss = %g%%, ratio %g\n', 100*(loss_total/init_loss_total), 100*(best_loss_total/init_loss_total), loss_total/best_loss_total);
    if(loss_total < best_loss_total)
        best_beta = beta;
        best_loss_total = loss_total;

        % update relative losses and plot the result
        epoch = epoch + 1;
        epoch_within_instance = epoch_within_instance + 1;
        for j=1:exp_count
            loss_acc{j}(epoch+1) = 100 *  (loss(j) / init_loss(j));
        end
        loss_acc_total(epoch+1) = 100 *  (loss_total / init_loss_total);
        % log the evolution of beta
        beta_history{optimization_instance}(epoch_within_instance+1,:) = beta;
        reactor_plot;
        drawnow;

        fprintf('Success! (epoch %d), Partial relative losses: [ ', epoch);
        for j=1:exp_count
            fprintf('%g%% ', loss_acc{j}(end));
        end
        fprintf('], total loss %g%%.\n', loss_acc_total(end));
    
    end

end


