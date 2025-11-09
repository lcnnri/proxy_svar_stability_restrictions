function print_console_summary(full, reg, stats_reg, dataSpec, cfg)
%PRINT_CONSOLE_SUMMARY
%   Print structural parameters (full sample and regimes)
%   and write the elasticity table to report/tables/table_spec5.txt.

    n_regimes = stats_reg.n_regimes;
    n = dataSpec.n;
    k = dataSpec.k;
    MBBcfg    = cfg.MBBcfg;
    modelSpec = cfg.modelSpec;
    alpha     = (1 - modelSpec.confidence)/2;
    qL = alpha*100; qU = (1-alpha)*100;
    alpha_est = alpha;                 % for readability
    n_bootstrap = MBBcfg.n_bootstrap;

    fprintf(' Display estimates...\n')
    fprintf('---------------------------------------------------------------\n');

    %% 1) Full-sample structural parameters
    fprintf('---------------------------------------------------------------\n');
    disp(' Printing structural estimates of Proxy-SVAR ')

    B_est = full.G(1:n,:);          % B (n x k)
    phi_est = full.G(n+1:end,:);    % phi (k x k)

    B_est_vec   = B_est(:);
    phi_est_vec = phi_est(:);
    B_est_B     = full.Boot_B;      % (n x k x B)
    phi_est_B   = full.Boot_phi;    % (k x k x B) if you saved it

    B_est_B_vec   = reshape(B_est_B, [], n_bootstrap);
    phi_est_B_vec = reshape(phi_est_B, [], n_bootstrap);

    B_diff   = B_est_B_vec   - B_est_vec;
    phi_diff = phi_est_B_vec - phi_est_vec;

    B_lb   = reshape(B_est_vec   - abs(prctile(B_diff',  qL))', n, k);
    B_ub   = reshape(B_est_vec   - (-abs(prctile(B_diff', qU)))', n, k);
    phi_lb = reshape(phi_est_vec - abs(prctile(phi_diff',  qL))', k, k);
    phi_ub = reshape(phi_est_vec - (-abs(prctile(phi_diff', qU)))', k, k);

    disp('Parameter estimates ==================');
    fprintf('\nPROXY SVAR (FULL-SAMPLE)===========================\n');
    for j=1:k
        fprintf('Column %d----\n',j);
        for i=1:n
            fprintf('B_%d%d:   &  %2.3f & (%2.3f, %2.3f)\n', ...
                i, j, B_est(i,j), B_lb(i,j), B_ub(i,j));
        end
        for i=1:k
            fprintf('phi_%d%d: &  %2.3f & (%2.3f, %2.3f)\n', ...
                i, j, phi_est(i,j), phi_lb(i,j), phi_ub(i,j));
        end
    end
    disp('======================================');

    %% 2) Regime-dependent structural parameters
    disp('Parameter estimates ===========================');
    fprintf('\nPROXY SVAR (STABILITY RESTRICTIONS)===========================\n');
    for r = 1:n_regimes
        fprintf('Regime %d ===========================\n', r);
        for i=1:n
            for j=1:n
                fprintf('G_%d%d: \\underset{(%2.3f , %2.3f)}{%2.3f}\n', ...
                    i,j, reg{r}.G_lb(i,j), reg{r}.G_ub(i,j), reg{r}.G(i,j));
            end
        end
        for i=1:k
            for j=1:n
                fprintf('phi_%d%d: \\underset{(%2.3f , %2.3f)}{%2.3f}\n', ...
                    i,j, reg{r}.G_lb(n+i,j), reg{r}.G_ub(n+i,j), reg{r}.G(n+i,j));
            end
        end
        for i=1:k
            for j=1:k
                fprintf('sigma_%d%d: \\underset{(%2.3f , %2.3f)}{%2.3f}\n', ...
                    i,j, reg{r}.G_lb(n+i,n+j), reg{r}.G_ub(n+i,n+j), reg{r}.G(n+i,n+j));
            end
        end
    end

    %% 3) Elasticities and bands (full sample and regimes)

    % ---------- Regime-specific (your formulas) ----------
    for ii = 1:n_regimes
        inv_G             = inv(reg{ii}.G);
        reg{ii}.phiG      = -inv_G(2,3)/inv_G(2,2);
        reg{ii}.phiTax    = -inv_G(1,3)/inv_G(1,1);
        reg{ii}.maxG      = max(reg{ii}.G_Multiplier);
        reg{ii}.maxTax    = max(reg{ii}.Tax_Multiplier);
        reg{ii}.lagmaxG   = find(reg{ii}.G_Multiplier==reg{ii}.maxG)-1;
        reg{ii}.lagmaxTax = find(reg{ii}.Tax_Multiplier==reg{ii}.maxTax)-1;

        reg{ii}.corrzT    = reg{ii}.G(4,1)/sqrt(sum(reg{ii}.G(4,:).^2));
        reg{ii}.corrzG    = reg{ii}.G(5,2)/sqrt(sum(reg{ii}.G(5,:).^2));
        reg{ii}.exozT     = reg{ii}.G(4,3)/sqrt(sum(reg{ii}.G(4,:).^2));
        reg{ii}.exozG     = reg{ii}.G(5,3)/sqrt(sum(reg{ii}.G(5,:).^2));

        % Bootstrap-based bands
        reg{ii}.phiG_boot      = NaN(1,n_bootstrap);
        reg{ii}.phiTax_boot    = NaN(1,n_bootstrap);
        reg{ii}.corrzT_boot    = NaN(1,n_bootstrap);
        reg{ii}.corrzG_boot    = NaN(1,n_bootstrap);
        reg{ii}.exozT_boot     = NaN(1,n_bootstrap);
        reg{ii}.exozG_boot     = NaN(1,n_bootstrap);

        for b = 1:n_bootstrap
            inv_G_boot               = inv(reg{ii}.G_boot(:,:,b));
            reg{ii}.phiG_boot(b)     = -inv_G_boot(2,3)/inv_G_boot(2,2);
            reg{ii}.phiTax_boot(b)   = -inv_G_boot(1,3)/inv_G_boot(1,1);

            reg{ii}.corrzT_boot(b)   = reg{ii}.G_boot(4,1,b) / sqrt(sum(reg{ii}.G_boot(4,:,b).^2));
            reg{ii}.corrzG_boot(b)   = reg{ii}.G_boot(5,2,b) / sqrt(sum(reg{ii}.G_boot(5,:,b).^2));
            reg{ii}.exozT_boot(b)    = reg{ii}.G_boot(4,3,b) / sqrt(sum(reg{ii}.G_boot(4,:,b).^2));
            reg{ii}.exozG_boot(b)    = reg{ii}.G_boot(5,3,b) / sqrt(sum(reg{ii}.G_boot(5,:,b).^2));
        end

        reg{ii}.phiG_lb   = reg{ii}.phiG ...
            - abs(prctile(reg{ii}.phiG_boot - reg{ii}.phiG, alpha_est*100));
        reg{ii}.phiG_ub   = reg{ii}.phiG ...
            - (-abs(prctile(reg{ii}.phiG_boot - reg{ii}.phiG, (1-alpha_est)*100)));

        reg{ii}.phiTax_lb = reg{ii}.phiTax ...
            - abs(prctile(reg{ii}.phiTax_boot - reg{ii}.phiTax, alpha_est*100));
        reg{ii}.phiTax_ub = reg{ii}.phiTax ...
            - (-abs(prctile(reg{ii}.phiTax_boot - reg{ii}.phiTax, (1-alpha_est)*100)));

        % Peaks: bands from IRF bands at peak horizon
        reg{ii}.maxG_lb   = reg{ii}.G_Multiplier_lb(reg{ii}.lagmaxG+1);
        reg{ii}.maxG_ub   = reg{ii}.G_Multiplier_ub(reg{ii}.lagmaxG+1);
        reg{ii}.maxTax_lb = reg{ii}.Tax_Multiplier_lb(reg{ii}.lagmaxTax+1);
        reg{ii}.maxTax_ub = reg{ii}.Tax_Multiplier_ub(reg{ii}.lagmaxTax+1);

        corrzT_diff = reg{ii}.corrzT_boot - reg{ii}.corrzT;
        reg{ii}.corrzT_lb = reg{ii}.corrzT ...
            - abs(prctile(corrzT_diff, alpha_est*100));
        reg{ii}.corrzT_ub = reg{ii}.corrzT ...
            - (-abs(prctile(corrzT_diff, (1-alpha_est)*100)));

        corrzG_diff = reg{ii}.corrzG_boot - reg{ii}.corrzG;
        reg{ii}.corrzG_lb = reg{ii}.corrzG ...
            - abs(prctile(corrzG_diff, alpha_est*100));
        reg{ii}.corrzG_ub = reg{ii}.corrzG ...
            - (-abs(prctile(corrzG_diff, (1-alpha_est)*100)));

        exozT_diff = reg{ii}.exozT_boot - reg{ii}.exozT;
        reg{ii}.exozT_lb  = reg{ii}.exozT ...
            - abs(prctile(exozT_diff, alpha_est*100));
        reg{ii}.exozT_ub  = reg{ii}.exozT ...
            - (-abs(prctile(exozT_diff, (1-alpha_est)*100)));

        exozG_diff = reg{ii}.exozG_boot - reg{ii}.exozG;
        reg{ii}.exozG_lb  = reg{ii}.exozG ...
            - abs(prctile(exozG_diff, alpha_est*100));
        reg{ii}.exozG_ub  = reg{ii}.exozG ...
            - (-abs(prctile(exozG_diff, (1-alpha_est)*100)));
    end

    % ---------- Full-sample elasticities ----------
    % A1 = (Sigma_yy^{-1} B)' as in your code
    A1        = (full.Sigma_Eta(1:n,1:n)\B_est)';  
    full_phiG   = -A1(2,3)/A1(2,2);
    full_phiTax = -A1(1,3)/A1(1,1);
    full_maxG   = max(full.G_Multiplier);
    full_maxTax = max(full.Tax_Multiplier);
    full_lagmaxG   = find(full.G_Multiplier == full_maxG)-1;
    full_lagmaxTax = find(full.Tax_Multiplier== full_maxTax)-1;

    % Corr(z, shock) with fiscal_shock = eps_y A1'
    fiscal_shock      = full.Errors(:,1:n)*A1';
    full_corr_ztax    = corr(fiscal_shock(:,1), full.Errors(:,4));
    full_corr_zg      = corr(fiscal_shock(:,2), full.Errors(:,5));

    % Bootstrap for full sample: phi and corr
    full_phiG_boot   = NaN(1,n_bootstrap);
    full_phiTax_boot = NaN(1,n_bootstrap);
    full_corr_ztax_boot = NaN(n_bootstrap,1);
    full_corr_zg_boot   = NaN(n_bootstrap,1);

    for b = 1:n_bootstrap
        B_b  = B_est_B(:,:,b);
        A1_b = (full.Sigma_Eta_boot(1:n,1:n,b)\B_b)';

        full_phiG_boot(b)   = -A1_b(2,3)/A1_b(2,2);
        full_phiTax_boot(b) = -A1_b(1,3)/A1_b(1,1);

        fiscal_shock_b = full.error(:,1:n,b)*A1_b';
        full_corr_ztax_boot(b) = corr(fiscal_shock_b(:,1), full.error(:,4,b));
        full_corr_zg_boot(b)   = corr(fiscal_shock_b(:,2), full.error(:,5,b));
    end

    % Bands
    full_phiG_lb   = full_phiG ...
        - abs(prctile(full_phiG_boot - full_phiG, alpha_est*100));
    full_phiG_ub   = full_phiG ...
        - (-abs(prctile(full_phiG_boot - full_phiG, (1-alpha_est)*100)));

    full_phiTax_lb = full_phiTax ...
        - abs(prctile(full_phiTax_boot - full_phiTax, alpha_est*100));
    full_phiTax_ub = full_phiTax ...
        - (-abs(prctile(full_phiTax_boot - full_phiTax, (1-alpha_est)*100)));

    full_maxG_lb   = full.G_Multiplier_lb(full_lagmaxG+1);
    full_maxG_ub   = full.G_Multiplier_ub(full_lagmaxG+1);
    full_maxTax_lb = full.Tax_Multiplier_lb(full_lagmaxTax+1);
    full_maxTax_ub = full.Tax_Multiplier_ub(full_lagmaxTax+1);

    qT_diff = full_corr_ztax_boot - full_corr_ztax;
    full_corr_ztax_lb = full_corr_ztax ...
        - abs(prctile(qT_diff, alpha_est*100));
    full_corr_ztax_ub = full_corr_ztax ...
        - (-abs(prctile(qT_diff, (1-alpha_est)*100)));

    qG_diff = full_corr_zg_boot - full_corr_zg;
    full_corr_zg_lb = full_corr_zg ...
        - abs(prctile(qG_diff, alpha_est*100));
    full_corr_zg_ub = full_corr_zg ...
        - (-abs(prctile(qG_diff, (1-alpha_est)*100)));

    % No exogeneity elasticities for full sample in your formulas -> set zero
    full_exozT = 0; full_exozT_lb = 0; full_exozT_ub = 0;
    full_exozG = 0; full_exozG_lb = 0; full_exozG_ub = 0;

    %% 4) Elasticity table (LaTeX, full + 2 regimes)
    fprintf('\nWriting elasticity table (full + 2 regimes)...\n');

    if ~isfolder('results/tables')
        mkdir('results/tables');
    end
    fid = fopen(fullfile('results','tables',['id' cfg.serial_id '_table_3.txt'] ),'w');

    desc = { ...
        '$\phi_{tax}$'; ...
        '$\mathcal{M}_{tax}$'; ...
        'relevance$_{tax}$ (\%)'; ...
        'contamination$_{tax}$ (\%)'; ...
        '$\phi_{g}$'; ...
        '$\mathcal{M}_{g}$'; ...
        'relevance$_{g}$ (\%)'; ...
        'contamination$_{g}$ (\%)'};

    us1 = '$\\underset{(%2.3f, %2.3f)}{%2.3f}$';
    us2 = '$\\underset{(%2.3f, %2.3f)}{%2.3f(%d)}$';

    % FULL / Regime 1 / Regime 2
    phiTax   = [full_phiTax;     reg{1}.phiTax;    reg{2}.phiTax];
    phiTax_l = [full_phiTax_lb;  reg{1}.phiTax_lb; reg{2}.phiTax_lb];
    phiTax_u = [full_phiTax_ub;  reg{1}.phiTax_ub; reg{2}.phiTax_ub];
    phiTaxtab = compose(us1, phiTax_l, phiTax_u, phiTax);

    maxTax   = [full_maxTax;     reg{1}.maxTax;    reg{2}.maxTax];
    maxTax_l = [full_maxTax_lb;  reg{1}.maxTax_lb; reg{2}.maxTax_lb];
    maxTax_u = [full_maxTax_ub;  reg{1}.maxTax_ub; reg{2}.maxTax_ub];
    maxTax_lag = [full_lagmaxTax; reg{1}.lagmaxTax; reg{2}.lagmaxTax];
    maxTaxtab  = compose(us2, maxTax_l, maxTax_u, maxTax, maxTax_lag);

    corrzT   = [full_corr_ztax;     reg{1}.corrzT;    reg{2}.corrzT];
    corrzT_l = [full_corr_ztax_lb;  reg{1}.corrzT_lb; reg{2}.corrzT_lb];
    corrzT_u = [full_corr_ztax_ub;  reg{1}.corrzT_ub; reg{2}.corrzT_ub];
    corrzTtab = compose(us1, corrzT_l, corrzT_u, corrzT);

    exoT   = [full_exozT;     reg{1}.exozT;    reg{2}.exozT];
    exoT_l = [full_exozT_lb;  reg{1}.exozT_lb; reg{2}.exozT_lb];
    exoT_u = [full_exozT_ub;  reg{1}.exozT_ub; reg{2}.exozT_ub];
    exoTtab = compose(us1, exoT_l, exoT_u, exoT);

    phiG   = [full_phiG;     reg{1}.phiG;     reg{2}.phiG];
    phiG_l = [full_phiG_lb;  reg{1}.phiG_lb;  reg{2}.phiG_lb];
    phiG_u = [full_phiG_ub;  reg{1}.phiG_ub;  reg{2}.phiG_ub];
    phiGtab = compose(us1, phiG_l, phiG_u, phiG);

    maxG   = [full_maxG;     reg{1}.maxG;     reg{2}.maxG];
    maxG_l = [full_maxG_lb;  reg{1}.maxG_lb;  reg{2}.maxG_lb];
    maxG_u = [full_maxG_ub;  reg{1}.maxG_ub;  reg{2}.maxG_ub];
    maxG_lag = [full_lagmaxG; reg{1}.lagmaxG; reg{2}.lagmaxG];
    maxGtab  = compose(us2, maxG_l, maxG_u, maxG, maxG_lag);

    corrzG   = [full_corr_zg;     reg{1}.corrzG;    reg{2}.corrzG];
    corrzG_l = [full_corr_zg_lb;  reg{1}.corrzG_lb; reg{2}.corrzG_lb];
    corrzG_u = [full_corr_zg_ub;  reg{1}.corrzG_ub; reg{2}.corrzG_ub];
    corrzGtab = compose(us1, corrzG_l, corrzG_u, corrzG);

    exoG   = [full_exozG;     reg{1}.exozG;    reg{2}.exozG];
    exoG_l = [full_exozG_lb;  reg{1}.exozG_lb; reg{2}.exozG_lb];
    exoG_u = [full_exozG_ub;  reg{1}.exozG_ub; reg{2}.exozG_ub];
    exoGtab = compose(us1, exoG_l, exoG_u, exoG);

    % helper for one row
    function print_row(lbl, cells)
        fprintf(fid, '%s & ', lbl);
        for jj = 1:numel(cells)
            if jj < numel(cells)
                fprintf(fid, '%s & ', cells{jj});
            else
                fprintf(fid, '%s ', cells{jj});
            end
        end
        fprintf(fid, '%s\n','\\');
    end

    print_row(desc{1}, phiTaxtab);
    print_row(desc{2}, maxTaxtab);
    print_row(desc{3}, corrzTtab);
    print_row(desc{4}, exoTtab);
    print_row(desc{5}, phiGtab);
    print_row(desc{6}, maxGtab);
    print_row(desc{7}, corrzGtab);
    print_row(desc{8}, exoGtab);

    fclose(fid);
    fprintf('Elasticity table written to results/tables/%s\n',['id' cfg.serial_id '_table_3.txt']);
end

