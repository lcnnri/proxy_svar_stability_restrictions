function save_results(full, reg, stats_reg, dataSpec, cfg)

if ~isfolder('results/workspace')
    mkdir('results/workspace');
end

save(fullfile('results','workspace', ...
     sprintf('id%s_fiscal_regimes.mat', cfg.serial_id)), ...
     'full','reg','stats_reg','dataSpec','cfg', '-v7.3');

end