function [full] = run_fullsample_proxysvar(dataSpec, cfg, rootDir)

% Unpack config
MBBcfg    = cfg.MBBcfg;
modelSpec = cfg.modelSpec;
alpha     = (1 - modelSpec.confidence)/2;
clock_    = cfg.clock_seed;
opt       = cfg.opt;



% Tax/GDP and G/GDP scaling constants (for multipliers)
scalingConstFullSample = struct();
scalingConstFullSample.TRY = 0.1823;
scalingConstFullSample.GY  = 0.2050;

fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' PROXY-SVAR ')
fprintf('---------------------------------------------------------------\n');
fprintf('\n');

%% 1) Reduced-form augmented constrained VAR (full sample) ================
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' 1) Estimate the reduced-form VAR...')

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
fprintf('\n');

[EstMdl, EstSE, logL, Errors] = estimate(Mdl, W);
% organize results
A_Const = EstMdl.Constant;
A = cell(1,p);
for L = 1:p, A{L} = EstMdl.AR{L}; end
Sigma_Eta = EstMdl.Covariance;

% First-stage F-statistics
est_lm = fitlm(W(p+1:end, n+1), Errors(:,1), 'Intercept', false);
[EstCoeffCov,se1,coeff1] = hac(W(p+1:end, n+1), Errors(:,1), 'Intercept', false, 'Display','off'); 
F_homo = (est_lm.Coefficients.Estimate/est_lm.Coefficients.SE)^2;
F_hete   = (coeff1/se1)^2;
fprintf('First-stage F (homoskedastic): %.2f | F (heteroskedasticy-robust): %.2f\n', F_homo, F_hete);

% Store reduced form objects
full = struct();
full.W = W;
full.A_Const = A_Const;
full.A = A;
full.Sigma_Eta = Sigma_Eta;
full.Errors = Errors;
full.logL = logL;
full.F_homo = F_homo;
full.F_hete = F_hete;

% Companion for IRFs
AL = [];
for L = 1:p, AL = [AL, A{L}]; end % concatenate coefficient matrices
Sj = [eye(M), zeros(M, M*(p-1))];
full.comp_matr = [AL; eye(M*(p-1)), zeros(M*(p-1),M)];

%% 2) Bootstrap the reduced-form ==========================================
% Perform bootstrap of the reduced-form augmented constrained VAR. If
% MBBcfg.MBB=true, use the specified size block; else set the size_block to
% 1. The latter choice collapses to the iid bootstrap. 
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' 2)Bootstrap the reduced-form VAR (MBB)... ')

if MBBcfg.MBB
    size_block = MBBcfg.size_block;
else
    size_block = 1; % iid bootstrap
end

bootFile = fullfile(rootDir,'FullSample_ProxySVAR_boot.mat');
if cfg.use_boot_cache && isfile(bootFile)
    S = load(bootFile); boot_Mdl = S.boot_Mdl;
else
    if MBBcfg.verbose, disp('Bootstrapping reduced-form VAR (MBB)...'); end
    [boot_Mdl, V_MBB] = do_VAR_MBB(full.Errors, size_block, full.W, full.A, full.A_Const, Mdl, ...
                                   MBBcfg.n_bootstrap, MBBcfg.verbose, MBBcfg.use_parallel, MBBcfg.num_workers,clock_);
    boot_Mdl.V_MBB = V_MBB;
    save(bootFile,'boot_Mdl');
    disp('Bootstrapping reduced-form VAR finished!');
end

% Store boot stuff
full.Sigma_Eta_boot = boot_Mdl.Sigma_eta;
full.Boot_A         = boot_Mdl.A;
full.Boot_Const     = boot_Mdl.Const;
full.comp_matr_boot = boot_Mdl.comp_matr;
full.error          = boot_Mdl.error;
full.V_MBB          = boot_Mdl.V_MBB;


%% 3) Proxy-SVAR estimation via minimum distance (point-identified) =======
% Proxy-SVAR, 2 instruments for 2 shocks
% Estimation on the full sample follows Angelini and Fanelli (2019, JAE)

% Preparing for the proxy-SVAR estimation... 
%   vR      (M x k):    mask containing equaility restrictions
%                       NaNs indicate free parameters (to estimate)
%   theta0 (Mk-nres x 1):  a vector of starting values for the free parameters
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' 3) Estimation of structural parameters of the Proxy-SVAR... ')


vR = [ nan nan;          % row 1..M, col 1..k (k instruments)
       nan nan;
       nan nan;
       nan nan;
         0 nan];         % zero restriction

vg0 = [ 0.0120  0.0018;  % starting values, shape (M x k)
        0.0055  0.0124;
       -0.0095  0.0025;
        0.0428  0.0000;
        0.0000  0.0129];
vg0 = vg0(1:M, 1:k); % unused

theta0 = vR;
theta0(isnan(vR)) = vg0(isnan(vR));
theta0 = theta0(isnan(vR(:)));      % free params vectorized

% Objective wrapper
obj_MD = @(theta,Sigma) md_partial(theta, Sigma, vR, k, k, 'V_MBB', full.V_MBB);

% Estimate on sample covariance
[Param, fval_MD] = fminunc(@(theta)obj_MD(theta, full.Sigma_Eta), theta0); 

% Recompose G = [B(:,1:k); phi(:,1:k)]
G_est = reshape(vR, M, k);
G_est(isnan(vR)) = Param;
B_est   = G_est(1:n,:);     % structural parameters
phi_est = G_est(n+1:end,:); % relevance of instruments


% Self-sign normalization
for i=1:k
    if B_est(i,i) < 0
        B_est(:,i)   = -B_est(:,i);
        phi_est(:,i) = -phi_est(:,i);
    end
end

% store point-identified part of impact matrix
full.G = [B_est; phi_est];


if size(phi_est,1)>size(B_est,2) % if there is more instruments than identified shocks
    % Overidentification test
    [fval, AVAR] = md_partial(Param, full.Sigma_Eta, vR, k, k, 'V_MBB', full.V_MBB);
    full.se = sqrt(diag(AVAR)/T);
    JVal    = T * fval;
    df_J    = numel(vR)-sum(isnan(vR(:)));
    if df_J > 0
        pJ = 1-chi2cdf(JVal, df_J);
        fprintf('Full-Sample Over-ID J: %.3f  p=%.3f  df=%d\n', JVal, pJ, df_J);
    end
end

%% 4) Bootstrap structural parameters and IRFs ============================
fprintf('\n');
fprintf('---------------------------------------------------------------\n');
disp(' 4) Bootstrap the structural parameters of the Proxy-SVAR... ')
if MBBcfg.verbose, disp('Bootstrapping structural parameters and IRFs (Full Sample)...'); end
[Boot_Param, Boot_Distance] = bootstrap_SVAR( ... 
    MBBcfg.n_bootstrap, @(ths,Sg)obj_MD(ths,Sg), Param, full.Sigma_Eta_boot, ...
    MBBcfg.verbose, MBBcfg.use_parallel, MBBcfg.num_workers,clock_);
if MBBcfg.verbose, disp('Bootstrapping of structural parameters is finished.'); end


B_est_B  = zeros(n,k,MBBcfg.n_bootstrap);
phi_est_B= zeros(k,k,MBBcfg.n_bootstrap);
full.se_bs = zeros(numel(Param), MBBcfg.n_bootstrap);

for b=1:MBBcfg.n_bootstrap
    [~, Avar_b] = obj_MD(Boot_Param(:,b), full.Sigma_Eta_boot(:,:,b));
    full.se_bs(:,b) = sqrt(diag(Avar_b)/T);

    Gb = reshape(vR, M, k);
    Gb(isnan(vR)) = Boot_Param(:,b);
    Bb   = Gb(1:n,:); 
    phib = Gb(n+1:end,:);

    % sign normalization
    for i=1:k
        if Bb(i,i) < 0
            Bb(:,i)   = -Bb(:,i);
            phib(:,i) = -phib(:,i);
        end
    end

    B_est_B(:,:,b)   = Bb;
    phi_est_B(:,:,b) = phib;
end

full.Boot_B = B_est_B;
full.Boot_phi = phi_est_B;

Multivariate_b_paper = [squeeze(B_est_B(1,1,1:15))';
                        squeeze(B_est_B(2,1,1:15))';
                        squeeze(B_est_B(3,1,1:15))';
                        squeeze(B_est_B(1,2,1:15))';
                        squeeze(B_est_B(2,2,1:15))';
                        squeeze(B_est_B(3,2,1:15))';
                            ];

fprintf(['Run bootstrap diagnostic test of identification\n ' ...
    'Refs: Angelini, Cavaliere, and Fanelli (2025, JoE), and Cavaliere, Fanelli, and Georgiev (2025) \n']);
disp('The test is a Normality test on the bootstrap estimates');
[ ~ , full.DornikHansen_Multivariate_PAPER] = DorHanomunortest(Multivariate_b_paper'); 



fprintf('\n');
disp('Build the IRFs and their confidence intervals');
% Sample IRFs
IRF = zeros(M,k,horizons+1);
for h=0:horizons
    IRF(:,:,h+1) = Sj * (full.comp_matr^h) * Sj' * [B_est; phi_est];
end
full.IRF = IRF;

full.Tax_Multiplier = -(squeeze(IRF(3,1,:))) / scalingConstFullSample.TRY / B_est(1,1);
if k>1
    full.G_Multiplier  =  (squeeze(IRF(3,2,:))) / scalingConstFullSample.GY  / B_est(2,2);
end

% Bootstrap IRFs
IRF_B = zeros(M,k,horizons+1,MBBcfg.n_bootstrap);
Tax_Mult_B = zeros(horizons+1, MBBcfg.n_bootstrap);
if k>1, G_Mult_B = zeros(horizons+1, MBBcfg.n_bootstrap); end

for b=1:MBBcfg.n_bootstrap
    for h=0:horizons
        IRF_B(:,:,h+1,b) = Sj * (full.comp_matr_boot(:,:,b)^h) * Sj' * [B_est_B(:,:,b); phi_est_B(:,:,b)];
    end
    Tax_Mult_B(:,b) = -(squeeze(IRF_B(3,1,:,b))) / scalingConstFullSample.TRY / B_est_B(1,1,b);
    if k>1
        G_Mult_B(:,b) =  (squeeze(IRF_B(3,2,:,b))) / scalingConstFullSample.GY / B_est_B(2,2,b);
    end
end

% Percentile bands 
qL = alpha*100; qU = (1-alpha)*100;

dT = Tax_Mult_B - full.Tax_Multiplier;
full.Tax_Multiplier_lb = full.Tax_Multiplier - abs(prctile(dT', qL))';
full.Tax_Multiplier_ub = full.Tax_Multiplier - (-abs(prctile(dT', qU)))';

if k>1
    dG = G_Mult_B - full.G_Multiplier;
    full.G_Multiplier_lb = full.G_Multiplier - abs(prctile(dG', qL))';
    full.G_Multiplier_ub = full.G_Multiplier - (-abs(prctile(dG', qU)))';
end





end