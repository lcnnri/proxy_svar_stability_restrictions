function cfg = default_config_fiscal()

cfg = struct();
cfg.serial_id           = '1'; % for archive purposes
% multistart:
cfg.multistart          = true; % multistart objective
cfg.n_starting_points   = 50;   % number of starting points for multistart
cfg.parallel_multistart = true; % use parallel
% plots:
cfg.plot_irf            = true;

% bootstrap:
cfg.MBBcfg = struct();
cfg.MBBcfg.n_bootstrap  = 4999; % number of bootstrap replications
cfg.MBBcfg.use_parallel = true; % run bootstrap in parallel
cfg.MBBcfg.verbose      = true; % print progress bar during bootstrap
cfg.MBBcfg.num_workers  = 8;    % max number of parallel workers for bootstrap
cfg.MBBcfg.MBB          = true; % true MBB, else iid bootstrap
cfg.MBBcfg.size_block   = 4; % block size for MBB
% use previous bootstrapped reduced form
cfg.use_boot_cache      = false; % true: does not perform VAR_MBB, it loads 
%                                   reduced-form bootstraps of the VAR


% model specification
cfg.modelSpec = struct();
cfg.modelSpec.p          = 4;    % VAR lag order
cfg.modelSpec.confidence = 0.68; % level of confidence for IRFs CI must lie in (0,1)
cfg.modelSpec.horizons   = 20;   % IRFs horizon

% reproducibility
cfg.clock_seed = 333*2; % set seed 


% Optim options (used across fminunc calls)
cfg.opt = optimset('MaxFunEvals',2e5,'TolFun',1e-12,'MaxIter',2e5,'TolX',1e-10,'Display','off');

end
