function [loss, Vg, S, Avar,  J, loss_w_reg] = ...
    hetcmd(theta,Shat,  vRG, vRL,V_MBB, const_irf)
%HETCMD computes the Classical Minimum-Distance estimator
%               with unconditional heteroscedasticity.
% Computes the classical minimum distance approach for estimation of the
%   structural parameters of a square system with parameter changes using
%   the covariances. References: Angelini, Fanelli & Neri (2025),
%   Lanne & Lutkepohl (2008), Bachiocchi Fanelli (2015), Angelini & Fanelli (2019).
%
% The baseline system is given by u = G * e. Examples:
%       Constant loading:
%           reg1 := G * G'
%           reg2 := G * L1 * G'
%           reg3 := G * L2 * G'
%           etc ...
%       Time-varying irf:
%       The structure of the volatility regimes is for
%           reg1 := G * G'
%           reg2 := (G + L1) * (G + L1)'
%           reg3 := (G + L1 + L2) * (G + L1 + L2)'
%           etc ...
%
%   Inputs:
%       theta       (double)     vector of structural parameters of size n^2*nsig
%       Shat        (nxnxr)  three-dimensional matrix of reduced-form covariance
%                                   matrices. Shat(:,:,1) related to first regime,
%                                   Shat(:,:,2) to the second regime etc.
%                                   otherwise
%       vRG         (n^2x1)       vector of restricted elements in G
%                                   (baseline regime)
%       vRL         (n^2(r-1)x1)  vector of restricted elements in L (every
%                                   regime >1
%       const_irf   (logical)     constant irf true (default), time-varying
%
%   Outputs:
%       loss       (scalar)       value of objective function to minimize
%       Vg         (double)       inverse of weighting matrix
%       S          (double)       local bias sensitivity to estimation error
%       Avar       (double)       asymptotic covariance matrix of parameters
%       J          (double)       Jacobian matrix of moment conditions
%       loss_w_reg (double)       loss value within each regime
%
% usage:
%   from bivariate VAR with 2 regimes 
%   Given reduced-form covariances of VAR: Shat(:,:,1), Shat(:,:,2)
%   Given estimate of variance of blkdiag(Shat(:,:,1),Shat(:,:,2)): V_MBB
%
%   vRG = nan(4,1); vRL = diag(nan(2,1));
%   theta0 = zeros(4+2,1);
%   hetcmd(theta,Shat,  vRG, vRL,V_MBB, const_irf)
% References: 
% - G. Angelini, L. Fanelli, and L. Neri (2025) "Invalid proxies 
%   and volatility changes" 
% - E. Bacchiocchi and L. Fanelli (2015, OBES) 
% - M. Lanne and H. Lutkepohl
% LN 2023
%=========================================================================%
arguments
    theta (:,1) double
    Shat (:,:,:) double
    vRG (:,:) double    = nan(size(Shat,1)^2,1);
    vRL (:,:,:) double  = nan(size(Shat,1)^2*(size(Shat,3)-1),1);
    V_MBB (:,:) double  = eye(size(Shat,1)*size(Shat,3));
    const_irf (1,1) logical=false;
end
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


[n,~,nsig] = size(Shat);    % number of covariances
nb = nsig-1;                % number of breaks
npars = numel(theta);       % number of parameters

% -
vR = [vRG;vRL]; % collect all restrictions
assert(sum(isnan(vR))==npars, ...
    ['The number of free parameters in vRG and vRL does not coincide with ' ...
    '    the number of elements in theta'])
thetalong=vR;
thetalong(isnan(vR)) = theta; % unpack parameters to estimate
theta = thetalong; 
thetalong=[]; %free space (not a big deal)

% -
G   = reshape(theta(1:n^2),n,n);
Sig = zeros(n,n,nsig);
D   = full(dupl(n)); % duplication matrix
Dp  = pinv(D);       % elimination matrix

if const_irf % Lanne & Lutkepohl
    F   = reshape(theta(n^2+1:end),n,n,nb);
    Sig = compute_sigmas(G,F,const_irf,nb); % Sig = G*F*G'
elseif ~const_irf % Bacchiocchi & Fanelli
    L   = reshape(theta(n^2+1:end),n,n,nb);
    Sig = compute_sigmas(G,L,const_irf,nb); % Sig = (G+L)*(G+L)'
end

Pgtr = zeros(nsig, n*(n+1)/2); % permuted vector of moments
for b=1:nsig
    Pgtr(b,:) = vech(Shat(:,:,b) - Sig(:,:,b))';
end
g = reshape(Pgtr',(n*(n+1)/2)*nsig,1); % vector of moments of correct size
Pgtr=[];


% Weighting matrix calculation / avoiding inverting big matrices
if ~isempty(V_MBB) % alpha-mixing residuals
    iVsig = V_MBB^-1;
else % iid residuals
    iVsig  = (1/2) * D' * kron(eye(n)/Shat(:,:,1),eye(n)/Shat(:,:,1)) * D;
    for m=2:nsig
        iVsig  = blkdiag(iVsig, ...
            (1/2) * D' * kron(eye(n)/Shat(:,:,m), eye(n)/Shat(:,:,m)) * D);
    end
end
iVg = iVsig;% weighting matrix, Vg^(-1)

loss = g' * iVg * g; % loss function

% % 
if nargout>1
    Vsig  = 2 * Dp * kron(Shat(:,:,1),Shat(:,:,1)) * Dp';
    for m=2:nsig
        Vsig  = blkdiag(Vsig, ...
            2 * Dp * kron(Shat(:,:,m),Shat(:,:,m)) * Dp');
    end
    Vg   = Vsig; % covariance of g
end
% % 
if nargout>=3 % compute sensitivities and informativeness
    
    if const_irf
        lam = F;
    elseif ~const_irf
        lam = L;
    end
    J  = compute_jacobian(G,lam,const_irf,n,nb,vRG,vRL); % get derivatives 
    %                                                      of moments
    
    % sens  = (J'* iVg * J) \ J' * iVg;    % sensitivity to finite-sample 
    sens  = pinv(J'* iVg * J) * J' * iVg;    % sensitivity to finite-sample 
    %                                          OLS estimation error 
    %                                          (use pinv, non-consequential 
    %                                          here, as AVAR is not used 
    %                                          anywhere)
    S       = sens * diag(sqrt(diag(Vg)));   % rescaled sensitivity
    Avar    = sens * Vg * sens'; % asymptotic variance of structural 
    %                              parameters: (G'WG)^(-1)
end

if nargout==6 % loss per regime (used for over-id restrictions test robust 
    %           to regime breaks)
    [loss_w_reg, ~] = equal_block_losses(g, iVg, n, nsig); % 
end


end %eof * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

function [Sig]= compute_sigmas(G,F,const_irf,nb) % * * * * * * * * * * * *
% COMPUTE_SIGMAS computes covariances using the matrices of structural
% parameters
%   Inputs:
%       G           (double) nxn    matrix of parameters
%       F           (double) nxnxnb matrix of parameters (if const_irf, it
%                                       should be diag)
%       const_irf   (logical)       if true use Lanne & Lutkepohl else 
%                                       Bacchiocchi & Fanelli (2015, OBES)
%       nb          (double)        number of breaks
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
Sig(:,:,1) = G * G';

if const_irf % Lanne & Lutkepohl
    for m=1:nb
        Sig(:,:,m+1)=  G*F(:,:,m)*G';
    end
elseif ~const_irf % Bacchiocchi & Fanelli
    for m=1:nb
        Sig(:,:,m+1)= (G+sum(F(:,:,1:m),3)) * ...
            (G+sum(F(:,:,1:m),3))';
    end
end

end % eof  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

function [J] = compute_jacobian(G,L,const_irf,n,nb, vRG,vRL) % * * * * *
% COMPUTE_JACOBIAN function to compute the Jacobians of the moment
%   conditions wrt the structural paramters of the model
%       Inputs:
%           G           (double) nxn     matrix of structural parameters
%           L           (double) nxnxnb  second matrices of structural
%                                           parameters
%           const_irf   (logical)          if true Lanne & Lutkepohl
%           n           (double)        equations in every system
%           vRG         (double) n^2x1   Restrictions on vec G
%           vRL         (double)         Restrictions on the other
%                                           matrices
%       Outputs:
%           J           (double) Jacobian
%
% References: G. Angelini, L. Fanelli and L. Neri (2025) "Invalid proxies 
% and volatility changes", Proposition 2, Eq. (20).
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
a = sum(isnan(vRG)); % number of free parameters in G
idxSg = find(isnan(vRG));  idxSl = find(isnan(vRL));
Sg = zeros((n)^2,a);
for j=1:size(Sg,2)
    Sg(idxSg(j),j)=1;
end
RL = reshape(vRL,n,n,nb); % accommodate multiple breaks
for m=1:nb
    vRL = vec(RL(:,:,m));
    b = sum(isnan(vRL)); % number of free parameters in L
    Sltmp = zeros((n)^2,b);
    for j=1:b
        Sltmp(idxSl(j),j)=1;
    end
    if m==1;Sl=Sltmp;else;Sl=blkdiag(Sl,Sltmp);end;Sltmp=[];
end
D = full(dupl(n)); Dp = pinv(D); nsig=nb+1; %number of Sigmas
Jtmp = kron(G,eye(n)); % first block of the Jacobian
if const_irf
    for m=1:nb % multiple breaks
        Jtmp = [Jtmp, zeros(n^2);
                kron(G*L(:,:,m),eye(n)), ...
                zeros(n^2,n^2*(m-1)), ...
                kron(G,G)];
    end
    Sl = 0.5*Sl; % rescale the selection matrix

elseif ~const_irf
    for m=1:nb % multiple breaks
        Jtmp = [Jtmp, zeros(m*n^2,n^2);
            repmat(kron(G+sum(L(:,:,1:m),3),eye(n)),1,m+1)];
    end
end

J = 2*kron(eye(nsig),Dp) * Jtmp * blkdiag(Sg,Sl);

end % eof * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%     J = 2*kron(eye(2),Dp) * ...
%         [kron(G,eye(n)), zeros(n^2);
%         kron(G*L,eye(n)), kron(G,G)] * ...
%         blkdiag(Sg,0.5*Sl);
%     J = 2*kron(eye(2),Dp) * ...
%         [kron(G,eye(n)), zeros(n^2); kron(G+L,eye(n)), kron(G+L,eye(n))] * ...
%         blkdiag(Sg,Sl);


function [losses, total] = equal_block_losses(g, iVg, n, r)
%EQUAL_BLOCK_LOSSES helper for computing the over-identifying restrictions 
% test statistic robust to regime breaks. We need the fval (loss) for each
% regime.
%   Inputs:
% g is the vector of moment conditions, and iVg is the weighting matrix.
% All blocks have size m = n(n+1)/2.
% r = number of losses (blocks) to compute.
%   Outputs: 
% losses is a vector of losses per regime. total is the total loss over the
% full sample. It matches loss in the main function.
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    m = n*(n+1)/2;
    if nargin < 4
        r = floor(numel(g)/m);        % infer r if not given
    end
    assert(numel(g) >= r*m, 'g too short for r blocks of size m');
    assert(all(size(iVg) == numel(g)), 'iVg size mismatch');

    losses = zeros(r,1);
    for k = 1:r
        idx = (k-1)*m + (1:m);
        gk   = g(idx);
        iVgk = iVg(idx,idx);
        losses(k) = gk' * iVgk * gk;
    end
    total = sum(losses);
end % eof * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *