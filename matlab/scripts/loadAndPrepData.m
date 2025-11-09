function [dataSpec,cfg, dataset] = loadAndPrepData(rootDir,dataSpec,cfg)
%LOADANDPREPDATA loads and prepare the data according to the dataSpec
%structure
% Outputs: 
%   dataSpec
%   cfg
%   dataset
%=========================================================================%
cfg.data_folder = fullfile(rootDir,'data',dataSpec.dataloc);
dataFile = fullfile(cfg.data_folder,'DATASET_AG.mat');
assert(isfile(dataFile), 'Data file not found: %s', dataFile);
S = load(dataFile);              % loads DATASET struct (LABEL, TSERIES, etc.)
dataset = S.DATASET;

% Locate Y columns
ixY = nan(1, dataSpec.n);
for j = 1:dataSpec.n
    hit = find(strcmpi(dataset.LABEL, dataSpec.columnnames{j}), 1);
    assert(~isempty(hit), 'Column not found: %s', dataSpec.columnnames{j});
    ixY(j) = hit;
end

% Locate instrument columns
ixIV = nan(1, dataSpec.k);
for j = 1:dataSpec.k
    hit = find(strcmpi(dataset.LABEL, dataSpec.instrumentnames{j}), 1);
    assert(~isempty(hit), 'Instrument not found: %s', dataSpec.instrumentnames{j});
    ixIV(j) = hit;
end

% Extract series (drop first row if it contains dates/headers)
years  = dataset.TSERIES(2:end,1);
ydata  = dataset.TSERIES(2:end, ixY);
z      = dataset.TSERIES(2:end, ixIV);
T      = size(ydata,1);


% Demean instruments on nonzero entries (MR convention)
if dataSpec.demean_instrument
    for j = 1:size(z,2)
        nz = z(:,j) ~= 0;
        if any(nz)
            z(nz,j) = z(nz,j) - mean(z(nz,j));
        end
    end
end

% Optional linear detrending of Y
if dataSpec.detrend_linear
    Y = ydata;
    X = [ones(T,1), (1:T)'];
    beta = (X' * X) \ (X' * Y);
    ydata = Y - X * beta;
end


dataSpec.years      = years;
dataSpec.ydata      = ydata;
dataSpec.z          = z;
dataSpec.T          = T;


end