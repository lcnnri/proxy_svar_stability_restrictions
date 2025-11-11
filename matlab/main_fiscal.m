% ========================================================================
%  Main replication script
%  Paper: "Invalid Proxies and Volatility Changes"
%  Authors: G. Angelini, L. Fanelli, L. Neri
% ========================================================================
% This script replicates the empirical illustration of
% "Invalid Proxies and Volatility Changes" (2025),
% by G. Angelini, L. Fanelli, and L. Neri.
%
% Data specification:
%   - Collects all information about the data, including file locations
%     and the variables used in the analysis.
%
% Configuration settings and data loading:
%
%   - Collects configuration settings for estimation and inference, such as:
%
%       * VAR lag order
%
%       * Number of bootstrap replications
%
%       * IRF horizon
%
%       * Multistart options for the objective function
%
%       * Parallel options
%
%   - Loads and prepares the data.
%
% Run proxy-SVARs:
%
%   a) [full] = run_fullsample_proxysvar(dataSpec, cfg, rootDir);
%      Estimates IRFs from a proxy-SVAR with k instruments and k shocks
%      (code can be generalized to k_1 instruments and k_2 < k_1 shocks).
%      Identification and estimation follow Angelini and Fanelli (2019, JAE).
%      Performs inference and computes MBB confidence intervals for
%      parameters and IRFs.
%
%   b) [reg, stats_reg] = run_regimes_proxysvar(dataSpec, DATASET, full, cfg);
%      Estimates time-varying IRFs from a proxy-SVAR identified using
%      stability restrictions as in Angelini, Fanelli, and Neri (2025, JBES).
%      Minimizes the CMD objective (eq. (23) in the paper) computed by:
%          hetcmd(theta, Shat, rG(:), rD(:), bigV_MBB)
%      where theta is the parameter vector, Shat the reduced-form
%      covariance matrices across regimes, rG the restrictions on the
%      reference-regime impact matrix, rD the stability restrictions for
%      subsequent regimes, and bigV_MBB the bootstrap covariance of the
%      reduced-form VAR error covariance matrix. hetcmd returns the
%      Jacobian of the moment conditions, used to check local identification
%      via Proposition 2.
%
% Plots and estimates:
%
%   - Produces the fiscal multipliers in Figure 1.
%
%   - Displays the summary of the estimated impact matrices as in Table S.5.
%
%   - Saves a text file with Table 3 to be copied into LaTeX.
%
% Save results:
%
%   - Saves all results to results\
%
% The code has been run and tested with the following system:
%
% ========================================================================
%
% MATLAB Version: 25.2.0.3042426 (R2025b) Update 1
% Operating System: Microsoft Windows 11 Home Version 10.0 (Build 26100)
% Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM)
%               64-Bit Server VM mixed mode
%
% ========================================================================
%
% MATLAB                                             Version 25.2 (R2025b)
%
% Required toolboxes:
%
%   - Optimization Toolbox
%
%   - Statistics and Machine Learning Toolbox
%
%   - Global Optimization Toolbox
%
%   - Econometrics Toolbox
%
%==========================================================================




clearvars; close all; clc;
tic;





% --- locate directories ---
thisFile    = mfilename('fullpath');
matlabDir   = fileparts(thisFile);
scriptsDir  = fullfile(matlabDir,'scripts');
funcDir     = fullfile(matlabDir,'functions');
dataDir     = fullfile(matlabDir,'data');
resultsDir  = fullfile(matlabDir,'results');
logDir      = fullfile(resultsDir,'logs');

addpath(scriptsDir, funcDir);

if ~exist(logDir,'dir'); mkdir(logDir); end

% --- start diary ---
logFileName = ['console_log_' datestr(now,'yyyymmdd_HHMMSS') '.txt'];
logFile = fullfile(logDir,logFileName );
if exist(logFile,'file'), delete(logFile); end
diary(logFile);
diary on;



fprintf('\n==============================================================\n');
fprintf("'Invalid proxies and volatility changes'\n")
fprintf(' Starting replication at %s\n', datestr(now));
fprintf(' Log file: %s\n', logFileName);
fprintf('==============================================================\n\n');


matlab_version_check;

%% Data specifications
dataSpec = struct();
dataSpec.dataloc= 'Fiscal/Caldara Kamps';
dataSpec.columnnames      = {'TAX_S','G_S','GDP_S'};   % endogenous Y vars
dataSpec.instrumentnames  = {'TAXNARRATIVE','AG'};     % external instruments
dataSpec.time             = 'Quarter';
dataSpec.demean_instrument= true;                       % MR-style demeaning
dataSpec.detrend_linear   = false;                      % linear detrend Y
dataSpec.n                = numel(dataSpec.columnnames);
dataSpec.k                = numel(dataSpec.instrumentnames);
dataSpec.M                = dataSpec.n + dataSpec.k;

rootDir = pwd;

%% Configurations settings and load data 
cfg = default_config_fiscal();    

% Echo settings
fprintf('--- Settings --------------------------------------------------\n');
fprintf('- VAR lags: %d\n', cfg.modelSpec.p);
fprintf('- Confidence level: %.2f\n', cfg.modelSpec.confidence);
fprintf('- Variables (Y): %s\n', strjoin(dataSpec.columnnames, ', '));
fprintf('- Instruments (Z): %s\n', strjoin(dataSpec.instrumentnames, ', '));
fprintf('- Time unit: %s\n', dataSpec.time);
fprintf('- IRF horizons: %d\n', cfg.modelSpec.horizons);


[dataSpec, cfg, DATASET] = loadAndPrepData(rootDir, dataSpec, cfg); % load data

%% Run proxy-SVARs
[full] = run_fullsample_proxysvar(dataSpec, cfg, rootDir);      % full-sample proxy-SVAR
[reg, stats_reg]  = run_regimes_proxysvar(dataSpec, DATASET, ... 
                                          full, cfg);           % proxy-SVAR with stability restrictions 

%% Plots and estimates
plot_multipliers(full, reg, dataSpec, stats_reg, cfg);                       % plots
print_console_summary(full, reg, stats_reg,dataSpec,cfg);                      % print estimates to screen

%% Save results
save_results(full, reg, stats_reg, dataSpec, cfg);


% End of script
fprintf('\n==============================================================\n');
fprintf(' Replication finished at %s\n', datestr(now));
fprintf('==============================================================\n');

diary off;
toc;