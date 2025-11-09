function [reg, stats_reg] = run_regimes_proxysvar(dataSpec, DATASET, full, cfg)

% Unpack config
MBBcfg    = cfg.MBBcfg;
modelSpec = cfg.modelSpec;
alpha     = (1 - modelSpec.confidence)/2;
clock_    = cfg.clock_seed;
opt = cfg.opt;

fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' PROXY-SVAR WITH STABILITY RESTRICTIONS')
fprintf('---------------------------------------------------------------\n');
fprintf('\n');

fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' 1) Estimate the reduced-form VAR (+ MBB)...')

n = dataSpec.n; k = dataSpec.k; M = dataSpec.M; p = modelSpec.p;
T = dataSpec.T; 
horizons = modelSpec.horizons;
ydata = dataSpec.ydata; z = dataSpec.z; 

% Build masks for augmented constrained VAR:
%  W = [Y', Z']'
% The proxies are already shocks, hence mean zero
VAR_Const = [nan(n,1); zeros(k,1)];  % constant in Y-eqs only
Ik = eye(k); Ik(Ik==1) = nan;
A_mask = [nan(M,n), [zeros(n,k); Ik]]; % AR part + contemporaneous Z in Z-eqs
VAR_A = cell(1,p);
for L = 1:p, VAR_A{L} = A_mask; end

% Define model
Mdl = varm('Constant',VAR_Const,'AR',VAR_A);
W   = [ydata z];

% Estimate augmented-constrained VAR
% %  Fixed break date (uncomment) -----------------------------------------
% T_B = 1983.25; % 1983Q2
% fprintf('\nBreak date (fixed): %4.2f\n', T_B);
% tau = find(years == T_B);
% tau0 = tau; 
% -------------------------------------------------------------------------
disp(' 1.a) Estimate break date and reduced-form VARs...')
% a) Estimate break date on the VAR's endogenous block only (done in the paper)
Mdl_rest = varm('Constant', VAR_Const(1:n), ...
                'AR', cellfun(@(A) A(1:n,1:n), VAR_A, 'UniformOutput', false));
% Candidate window around beginning of 1980s: [109:152] 
tau_index = estimate_break_point(full.W(:,1:n), Mdl_rest, 109:152);
T_B = dataSpec.years(tau_index); 
tau = tau_index / size(full.W,1); 
fprintf('Estimated break date: %dQ%d (tau %2.3f)\n', floor(T_B), (T_B-floor(T_B))/0.25+1, tau);

n_regimes = 2;

ttau = [0; tau_index; size(ydata,1)];
t_B = diff(ttau); % obs per regime

reg = cell(1,n_regimes); % init regime-specific structures
for r = 1:n_regimes 
    reg{r} = struct();
    reg{r}.W = full.W(ttau(r)+1:ttau(r+1),:); % split the data into regimes
end

fprintf(' 1.b) Estimate VAR with %d regimes...\n', n_regimes)
fprintf('   and bootstrap\n')

%  Reduced-form per-regime + bootstrap
for r = 1:n_regimes
    fprintf('       Estimate the VAR for regime %d...\n',r)
    [reg{r}.EstMdl, reg{r}.EstSE, reg{r}.logL, reg{r}.Errors] = ...
        estimate(varm('Constant',VAR_Const,'AR',VAR_A), reg{r}.W);
    reg{r}.A_Const = reg{r}.EstMdl.Constant;
    reg{r}.A = cell(1,p);
    for L=1:p, reg{r}.A{L} = reg{r}.EstMdl.AR{L}; end
    reg{r}.Sigma_Eta = reg{r}.EstMdl.Covariance;

    % Companion form
    ALr = [];
    for L=1:p, ALr = [ALr, reg{r}.A{L}]; end 
    reg{r}.comp_matr = [ALr; eye(M*(p-1)), zeros(M*(p-1),M)];

    if MBBcfg.MBB
    size_block = MBBcfg.size_block;
    else
        size_block = 1; % iid bootstrap
    end


    fprintf('       Bootstrap the VAR for regime %d...\n',r)
    cacheFile = sprintf('reg%d_boot.mat', r);
    if cfg.use_boot_cache && isfile(cacheFile)
        S = load(cacheFile); bM = S.bM;
    else
        fprintf('Bootstrapping reduced-form VAR (MBB) regime %d...\n', r);
        % Run bootstrap of VAR Bootstrap for each regime
        [bM, reg{r}.V_MBB] = do_VAR_MBB(reg{r}.Errors, size_block, reg{r}.W, ...
                                          reg{r}.A, reg{r}.A_Const, ...
                                          varm('Constant',VAR_Const,'AR',VAR_A), ...
                                          MBBcfg.n_bootstrap, MBBcfg.verbose, ...
                                          MBBcfg.use_parallel, MBBcfg.num_workers,clock_);
        fprintf('Bootstrapping reduced-form VAR for regime %d finished!\n', r);
        bM.V_MBB = reg{r}.V_MBB; % covariance of red-form covariances
        save(cacheFile, 'bM','-v7.3');
    end
    if ~isfield(bM,'V_MBB'), bM.V_MBB = reg{r}.V_MBB; end

    reg{r}.V_MBB         = bM.V_MBB;
    reg{r}.Sigma_Eta_boot= bM.Sigma_eta;
    reg{r}.Boot_A        = bM.A;
    reg{r}.Boot_Const    = bM.Const;
    reg{r}.comp_matr_boot= bM.comp_matr;
    reg{r}.error_boot    = bM.error;
end

% Stack V_MBB. Used for stability-restrictions CMD objective
V_MBB_blk = blkdiag(reg{1}.V_MBB, reg{2}.V_MBB);

% e) Regime-specific scaling for multipliers
% Compute TRY and GY per regime from original DATASET
scalingConstReg.TRY = zeros(1,n_regimes); % init tax revenues
scalingConstReg.GY  = zeros(1,n_regimes); % init gov 
scalingConstReg.TRY(1) = 1/mean(exp(DATASET.TSERIES(ttau(1)+1:ttau(2),2) - DATASET.TSERIES(ttau(1)+1:ttau(2),3)));
scalingConstReg.GY(1)  = 1/mean(exp(DATASET.TSERIES(ttau(1)+1:ttau(2),2) - DATASET.TSERIES(ttau(1)+1:ttau(2),4)));
scalingConstReg.TRY(2) = 1/mean(exp(DATASET.TSERIES(ttau(2)+1:ttau(3),2) - DATASET.TSERIES(ttau(2)+1:ttau(3),3)));
scalingConstReg.GY(2)  = 1/mean(exp(DATASET.TSERIES(ttau(2)+1:ttau(3),2) - DATASET.TSERIES(ttau(2)+1:ttau(3),4)));

%% 2) Identification via stability restrictions
% Preparing for the proxy-SVAR estimation across regimes using
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' 2) Estimation of structural parameters of the Proxy-SVAR via stability restrictions... ')
%   rG      (M x M):    mask containing equaility restrictions on the
%                       first regime's impact matrix
%                       NaNs indicate free parameters (to estimate)
%   rD      (M x M):    mask containing stability restrictions on the
%                       second regime's Delta matrix
%                       NaNs indicate free parameters (to estimate)
%   theta0 ((MM-nres) x 1):    a vector of starting values for the free
%                       parameters

% rG and rD patterns (M x M). Keep your original structure.
rG = [nan nan nan 0   0;
      nan nan 0   0   0;
      nan nan nan 0   0;
      nan nan nan nan 0;
      0   nan nan nan nan];

rD = [nan nan nan 0   0;
      nan nan 0   0   0;
      0   nan nan 0   0;
      nan nan 0   0   0;
      0   nan 0   0   nan];

% Starting values on matrix G and Delta
G0 = [full.G, ...
      [0.0096  0       0;
       0       0       0;
       0.0091  0       0;
       0.0024  0.0794  0;
       0.0024  0       0.0041]];


D0 = [ 0.0032  -0.0022  -0.0139   0    0;
       0.0005  -0.0069   0        0    0;
       0       -0.0015  -0.0038   0    0;
       0.0481   0        0        0    0;
       0       -0.0067   0        0   -0.0025];

% Trim to MxM % unused
G0 = G0(1:M,1:M); 
D0 = D0(1:M,1:M);
rG = rG(1:M,1:M); 
rD = rD(1:M,1:M);

THETA0 = [G0, D0]; % all parameters
rTot   = [rG, rD]; % all restrictions

% Order condition
num_free_pars = sum(isnan(rG(:))) + sum(isnan(rD(:))); % # free parameters
num_equations = n_regimes * M*(M+1)/2; % # equations
fprintf('Order condition %s: equations=%d, free parameters=%d\n', ...
        ternary(num_free_pars>num_equations, 'NOT satisfied','satisfied'), ...
        num_equations, num_free_pars);

% Stack covariances
Shat = zeros(M,M,n_regimes);
for r=1:n_regimes, Shat(:,:,r) = reg{r}.Sigma_Eta; end

% Estimate hetero-CMD with/without multistart
theta0_vec = THETA0(isnan(rTot));
bigV_MBB=blkdiag(reg{1}.V_MBB, reg{2}.V_MBB);
if ~cfg.multistart
    [thetahat, fval_cmd, exitflag, info] = ...
        fminunc(@(x)hetcmd(x, Shat, rG(:), rD(:), bigV_MBB, false), ...
                theta0_vec, opt); 
else
    fprintf('Multistart estimation procedure with %d random starting points\n', cfg.n_starting_points);
    rng(clock_); % fix seed
    ms = MultiStart('UseParallel', cfg.parallel_multistart);
    problem = createOptimProblem('fminunc','objective', ...
        @(x)hetcmd(x, Shat,  rG(:), rD(:), bigV_MBB, false), ...
        'x0', theta0_vec, 'options', opt);
    [thetahat, fval_cmd, exitflag, info] = run(ms, problem, cfg.n_starting_points); 
end


%% 3) Identification heuristics, and over-id restrictions test
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
fprintf(' 3) Identification checks and over-id restrictions test...\n')

% Sandwich, losses, and diagnostics
[~, ~, S_cmd, Avar_cmd, J_cmd,losses] = hetcmd(thetahat, Shat, rG(:), rD(:), bigV_MBB, false);
loss_1 = losses(1); loss_2 = losses(2); % losses in each reagime
se_cmd = sqrt(diag(Avar_cmd) / size(full.W,1)); % asymptotic se

% identifiability checks
rcondnum = rcond(J_cmd'*J_cmd);
eigsJJ = eig(J_cmd'*J_cmd);
condnum = max(eigsJJ) / min(eigsJJ);
S_J = svd(J_cmd'*J_cmd);
minSV = min(S_J);

% overidentifying restrictions test 
if num_equations>num_free_pars
    Jtest_Val = t_B(1)*loss_1 + t_B(2)*loss_2;
    pvalJtest = 1 - chi2cdf(Jtest_Val, num_equations - num_free_pars);
    fprintf('Regime-Change:  Over-ID J (CMD): %.3f, p=%.3f\n', Jtest_Val, pvalJtest);
end
% Rank condition of Proposition 2
if rank(J_cmd) ~= num_free_pars
    error('Rank condition not satisfied: rank(J)=%d, needed=%d', rank(J_cmd), num_free_pars);
else 
    fprintf('Rank condition satisfied.\n');
end

% Recover G for each regime
rTot(isnan(rTot)) = thetahat;
G_est = rTot(1:M,1:M);
D_est = rTot(1:M,M+1:end);
reg{1}.G = G_est;
reg{2}.G = G_est + D_est;

% Sign normalization
for r=1:n_regimes
    for i=1:M
        if reg{r}.G(i,i) < 0
            reg{r}.G(:,i) = -reg{r}.G(:,i);
        end
    end
end

%% 4) Sample IRFs by regime
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
fprintf(' 4) Construct IRFs by regime...\n')
% Sample IRFs
Sj  = [eye(M), zeros(M, M*(p-1))];

for r=1:n_regimes
    reg{r}.IRF = zeros(M,M,horizons+1);
    for h=0:horizons
        reg{r}.IRF(:,:,h+1) = Sj * reg{r}.comp_matr^h * Sj' * reg{r}.G;
    end
    reg{r}.Tax_Multiplier = -(squeeze(reg{r}.IRF(3,1,:))) / scalingConstReg.TRY(r) / reg{r}.G(1,1);
    if k>1
        reg{r}.G_Multiplier =  (squeeze(reg{r}.IRF(3,2,:))) / scalingConstReg.GY(r)  / reg{r}.G(2,2);
    end
end


%% 5) Bootstrap of the structural parameters
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
fprintf(' 5) Bootstrap the structural parameters of the Proxy-SVAR...\n')
% Bootstrap IRFs
Sigma_Eta_Boot = zeros(M,M,MBBcfg.n_bootstrap,n_regimes);
for r=1:n_regimes
    Sigma_Eta_Boot(:,:,:,r) = reg{r}.Sigma_Eta_boot(:,:,1:MBBcfg.n_bootstrap);
end


se_bs = nan(numel(thetahat), MBBcfg.n_bootstrap);
rcondnum_boot = nan(MBBcfg.n_bootstrap,1);
mineig = nan(MBBcfg.n_bootstrap,1);



objHet_MD = @(ths,Sg)hetcmd(ths, Sg,  rG(:), rD(:), bigV_MBB, false); % CMD objective as a handle

if MBBcfg.verbose, disp('Bootstrapping structural parameters (with stability restrictions)...'); end
[Boot_Param, Boot_Distance] = bootstrap_SVAR( ...
    MBBcfg.n_bootstrap, @(ths,Sg)objHet_MD(ths,Sg), thetahat, Sigma_Eta_Boot, ...
    MBBcfg.verbose, MBBcfg.use_parallel, MBBcfg.num_workers,clock_,true);
if MBBcfg.verbose, disp('Bootstrapping of structural parameters is finished.'); end



for b=1:MBBcfg.n_bootstrap    
    th_b = Boot_Param(:,b);
    Shat_b = cat(3, reg{1}.Sigma_Eta_boot(:,:,b), reg{2}.Sigma_Eta_boot(:,:,b));
    
    % store min SVD for identifiability check
    [~, ~, ~, Avar_b,J_b] = hetcmd(th_b, Shat_b,  rG(:), rD(:),bigV_MBB); % get Jacobian
    se_bs(:,b) = sqrt(diag(Avar_b) / size(full.W,1)); % asymptotic ses
    S_J_b = svd(J_b'*J_b);
    mineig(b) = min(S_J_b);
    rcondnum_boot(b) = rcond(J_b'*J_b);

    rTot_b = [rG, rD];
    rTot_b(isnan(rTot_b)) = th_b;
    G_b = rTot_b(1:M,1:M);
    D_b = rTot_b(1:M,M+1:end);

    G1 = G_b; 
    G2 = G_b + D_b;
    % Sign normalization
    for i=1:M
        if G1(i,i) < 0, G1(:,i) = -G1(:,i); D_b(:,i) = -D_b(:,i); end
        if G2(i,i) < 0, G2(:,i) = -G2(:,i); end
    end
    G_boot{1}(:,:,b) = G1;
    G_boot{2}(:,:,b) = G2;
end
% G_boot{1} = G_boot1;
% G_boot{2} = G_boot2;
% clear G_boot1 G_boot2;

% Build IRFs for each bootstrap draw
for r=1:n_regimes
    reg{r}.G_boot = G_boot{r};
    reg{r}.IRF_boot = zeros(M,M,horizons+1,MBBcfg.n_bootstrap);
    reg{r}.Tax_Multiplier_boot = zeros(horizons+1, MBBcfg.n_bootstrap);
    if k>1, reg{r}.G_Multiplier_boot = zeros(horizons+1, MBBcfg.n_bootstrap); end

    for j=1:MBBcfg.n_bootstrap
        for h=0:horizons
            reg{r}.IRF_boot(:,:,h+1,j) = Sj * reg{r}.comp_matr_boot(:,:,j)^h * Sj' * reg{r}.G_boot(:,:,j);
        end
        reg{r}.Tax_Multiplier_boot(:,j) = -(squeeze(reg{r}.IRF_boot(3,1,:,j))) / scalingConstReg.TRY(r) / reg{r}.G_boot(1,1,j);
        if k>1
            reg{r}.G_Multiplier_boot(:,j) =  (squeeze(reg{r}.IRF_boot(3,2,:,j))) / scalingConstReg.GY(r)  / reg{r}.G_boot(2,2,j);
        end
    end
end

% Confidence intervals per regime
qL = alpha*100; qU = (1-alpha)*100;
for r=1:n_regimes
    reg{r}.G_Quantiles = reshape(reg{r}.G_boot, M*M, MBBcfg.n_bootstrap) - repmat(reg{r}.G(:),1,MBBcfg.n_bootstrap);
    reg{r}.G_lb = reshape(reg{r}.G(:) - abs(prctile(reg{r}.G_Quantiles', qL))', M, M);
    reg{r}.G_ub = reshape(reg{r}.G(:) - (-abs(prctile(reg{r}.G_Quantiles', qU)))', M, M);

    dTax = reg{r}.Tax_Multiplier_boot - reg{r}.Tax_Multiplier;
    reg{r}.Tax_Multiplier_lb = reg{r}.Tax_Multiplier - abs(prctile(dTax', qL))';
    reg{r}.Tax_Multiplier_ub = reg{r}.Tax_Multiplier - (-abs(prctile(dTax', qU)))';

    if k>1
        dG = reg{r}.G_Multiplier_boot - reg{r}.G_Multiplier;
        reg{r}.G_Multiplier_lb = reg{r}.G_Multiplier - abs(prctile(dG', qL))';
        reg{r}.G_Multiplier_ub = reg{r}.G_Multiplier - (-abs(prctile(dG', qU)))';
    end
end

stats_reg = struct();
stats_reg.n_regimes       = n_regimes;
stats_reg.mineig          = mineig;
stats_reg.rcondnum        = rcondnum;
stats_reg.rcondnum_boot   = rcondnum_boot;
stats_reg.se_bs           = se_bs;
stats_reg.t_B             = t_B;
stats_reg.T_B             = T_B;
stats_reg.rG              = rG;
stats_reg.rD              = rD;
stats_reg.G_est           = G_est;
stats_reg.D_est           = D_est;

end
