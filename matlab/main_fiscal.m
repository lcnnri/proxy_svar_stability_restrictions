%% ========================================================================
%  Main replication script
%  Paper: "Invalid Proxies and Volatility Changes"
%  Authors: G. Angelini, L. Fanelli, L. Neri
% ========================================================================

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