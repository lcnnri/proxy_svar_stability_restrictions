function cfg = default_config_fiscal()

cfg = struct();
cfg.serial_id           = '1';
cfg.multistart          = true;
cfg.n_starting_points   = 50;
cfg.parallel_multistart = true;
cfg.plot_irf            = true;
cfg.use_boot_cache      = true;

cfg.MBBcfg = struct();
cfg.MBBcfg.n_bootstrap  = 4999;
cfg.MBBcfg.use_parallel = true;
cfg.MBBcfg.verbose      = true;
cfg.MBBcfg.num_workers  = 8;
cfg.MBBcfg.MBB          = true;
cfg.MBBcfg.size_block   = 4;

cfg.modelSpec = struct();
cfg.modelSpec.p          = 4;
cfg.modelSpec.confidence = 0.68;
cfg.modelSpec.horizons   = 20;

cfg.clock_seed = 333*2;


% Optim options (used across fminunc calls)
cfg.opt = optimset('MaxFunEvals',2e5,'TolFun',1e-12,'MaxIter',2e5,'TolX',1e-10,'Display','off');

end
