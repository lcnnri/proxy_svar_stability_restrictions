function [khat] = estimate_break_point(W,Mdl,range)
%ESTIMATE_BREAK_POINT Function to estimate the break points in a VAR. 
%   Function that implements the procedure of Bai (2000).
% Inputs
%   W      : T-by-M matrix of observables (rows: time, cols: variables).
%   Mdl    : VAR model template (varm) with chosen lag order p and any
%            fixed options. Estimation is done via estimate(Mdl, subW).
%   varargin{1} (optional):
%            'range' = K-by-J integer matrix of candidate break indices.
%            - K = number of breaks under consideration.
%            - J = number of candidate configurations to evaluate.
%            Column j collects the K break indices (strictly increasing).
%            If omitted, the function defaults to all single-break
%            candidates at 1:T, which is not meaningful in practice; in
%            applications you should pass a trimmed grid (e.g. disallow
%            breaks too close to the sample ends or to each other).
%
% Output
%   khat   : K-by-1 vector of selected break indices (the column of 'range'
%            that maximizes the total log-likelihood). Empty if no valid
%            estimation was successful.

arguments
    W (:,:) double
    Mdl varm
    range (1,:) double = 1:length(W) % Fallback: evaluate single-break grid over 1:T.
end

warning('off') % % suppress estimation warnings (e.g. near-singular covariances likely in small samples/many lags)
T = length(W);

% Number of regimes implied by K breaks
% If 'range' is K-by-J, we have n_regimes = K+1
n_regimes = size(range,1) + 1;

% Preallocate storage for objective and chosen breaks across configurations
Lambda    = zeros(size(range,2),1);
k_store = zeros(size(range));

% Loop over candidate configurations (columns of 'range')
for i=1:numel(range(1,:))
    % Build the break vector k = [0; k_1; ...; k_K; T]
    k = [];
    for r=1:n_regimes-1
        k(r,1) = range(r,i);
    end
    k = [0; k; T];

    % Regime lengths T_k = k_{r+1} - k_r
    T_k  = k(2:end)-k(1:end-1); % number of observations in every regime

    % Split and estimate regime-by-regime
    for ii = 1:n_regimes % split sample
        % Subsample indices: (k(ii)+1) : k(ii+1)
        subsW = W(k(ii)+1:k(ii+1),:);
        err=false;
 
        % Try estimation of the VAR on this regime
        try
        [~, ~,logL,~] = ...
                estimate(Mdl,subsW);
        catch
            err = true;
        end
        
        % Record objective and breaks if all regimes were successfully estimated
        if ~err
            Lambda(i) = Lambda(i) + logL;
        end

    end

    k_store(:,i) = k(2:end-1);
end

% Select the configuration with the highest total log-likelihood
[~,idx] = max(Lambda);
khat    = k_store(:,idx);
end