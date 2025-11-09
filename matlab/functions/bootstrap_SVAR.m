function [b_theta_hat,b_loss] = bootstrap_SVAR(n_bootstrap,func,thetahat,Sigma,verbose,use_parallel,num_workers,rng_clock,opt_flag)
%BOOTSTRAP_SVAR This function bootstraps the structural parameters of a
%SVAR model, identified from its covariances, Sigma.
%   Detailed explanation goes here
arguments
    n_bootstrap (1,1) double
    func (:,:) function_handle
    thetahat (:,1) double
    Sigma (:,:,:,:) double
    verbose (1,1) logical = true
    use_parallel (1,1) logical = true
    num_workers (1,1) double = 4
    rng_clock (1,1) double = 0
    opt_flag (1,1) logical = false
end


if opt_flag
    options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-50,...
        'Display','off');
else
    options = optimset('Display','off');
end

b_theta_hat = zeros(numel(thetahat),n_bootstrap);
b_loss      = nan(n_bootstrap,1);
completed = 0;   % number of completed iterations

if ~use_parallel % not parallel
    rng(rng_clock);
    for b=1:n_bootstrap
        Sb = cat(3,Sigma(:,:,b,:)); 
        [b_theta_hat(:,b), b_loss(b)] = ...
            fminunc(@(theta)func(theta, Sb), thetahat, options);
        if verbose
            message;
        end
    end
elseif use_parallel % use parallel
    pool = gcp('nocreate');  % Get current parallel pool (if any), or [] if none
    if isempty(pool)
        pool = ... Start a new parallel pool with default settings
            parpool(num_workers);  
    end

    % If it is verbose, we print a progress bar
    if verbose
        PB = ProgressBar(n_bootstrap, taskname='SVAR MBB', ui='cli', no_log=false);
    end
    rng(rng_clock)
    parfor (b=1:n_bootstrap,num_workers)
        Sb = cat(3,Sigma(:,:,b,:)); 
        [b_theta_hat(:,b), b_loss(b)] = ...
            fminunc(@(theta)func(theta, Sb), thetahat, options);
        if verbose
            count(PB)
        end
    end
end
    function message(~)
        completed = completed + 1;
        if mod(completed, 50) == 0
            fprintf('Bootstrapping the SVAR,  %4.0f%s completed (%d iterations out of %d) \n', ...
                completed/n_bootstrap*100, '%',completed,n_bootstrap);
        end
    end

end
