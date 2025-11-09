function [boot_Mdl,V_MBB] = do_VAR_MBB(uhat,block_length,y,A,Const,Mdl,n_bootstrap,verbose,use_parallel,num_workers,rng_clock)
%DO_VAR_MBB Summary of this function goes here
%   Detailed explanation goes here
arguments
    uhat (:,:) double
    block_length (1,1) double 
    y (:,:) double
    A  (1,:) cell
    Const (:,1) double
    Mdl varm
    n_bootstrap (1,1) double
    verbose (1,1) logical
    use_parallel (1,1) logical
    num_workers (1,1) double = 4
    rng_clock (1,1) double = 0
end 

% checks
if size(uhat,1)<size(uhat,2); uhat=uhat'; end
if size(y,1)<size(y,2); y=y'; end
assert(length(uhat) == length(y)-numel(A), 'Length of Y should be T+p');


[T,n] = size(uhat);
p     = numel(A);

numblocks = T-block_length+1;
N         = ceil(T/block_length);

u_blocks = zeros(block_length, n, numblocks);
for i=1:numblocks
    u_blocks(:,:,i) = uhat(i:block_length+i-1,:);
end

u_tmp = zeros(block_length,n);
for i=1:block_length
    for j=1:n
        u_tmp(i,j)=mean(uhat(i:numblocks+i-1,j));
    end
end

numResample = N;
u_center = zeros(numResample*block_length,n);
for i=1:numResample
    u_center((i-1)*block_length+1:i*block_length,:) = u_tmp;
end

% initialize output matrices
boot_Const = zeros(n,1,n_bootstrap);
boot_A = zeros(n,n,n_bootstrap,p);
boot_Sigma_eta = zeros(n,n,n_bootstrap);
boot_Comp_Mat  = zeros(n*p,n*p,n_bootstrap);
boot_error = zeros(T-p,n,n_bootstrap);
boot_y = zeros(T,n,n_bootstrap);

ratio = N;
if rng_clock==0
    rng();
else
    rng(rng_clock)
end
index = randi(numblocks, n_bootstrap, ratio);

if ~use_parallel
    % start bootstrap
    for b=1:n_bootstrap
        if verbose
            if mod(b,50)==0
                fprintf('Bootstrapping reduced-form of the SVAR, sample n. % d\n', b);
            end
        end
        [b_Const, b_A, b_Sigma_eta, b_Comp_Mat, b_error] = ...
            MBB_iteration(u_blocks,u_center,T,n,p,y,Const,A,Mdl,block_length,numblocks,ratio,index(b,:));

        boot_Const(:,:,b)= b_Const;
        boot_A(:,:,b,:) = b_A;
        boot_Sigma_eta(:,:,b) = b_Sigma_eta;
        boot_Comp_Mat(:,:,b) = b_Comp_Mat;
        boot_error(:,:,b) = b_error;
    end
elseif use_parallel % parallel bootstrap
    if verbose
        PB = ProgressBar(n_bootstrap, taskname='VAR MBB', ui='cli', no_log=false);
    end
    parfor (b=1:n_bootstrap, num_workers) 
        [b_Const, b_A, b_Sigma_eta, b_Comp_Mat, b_error, b_y] = ...
            MBB_iteration(u_blocks,u_center,T,n,p,y,Const,A,Mdl,block_length,numblocks,ratio,index(b,:));

        boot_Const(:,:,b)= b_Const;
        boot_A(:,:,b,:) = b_A;
        boot_Sigma_eta(:,:,b) = b_Sigma_eta;
        boot_Comp_Mat(:,:,b) = b_Comp_Mat;
        boot_error(:,:,b) = b_error;
        boot_y(:,:,b) = b_y;
        if verbose
            count(PB);
        end
    end
end

for b=1:n_bootstrap
    for l=1:p
        tmpb_A{l}(:,:,b) = boot_A(:,:,b,l);
    end
end
boot_A = tmpb_A;

% store
boot_Mdl.numblocks  = numblocks;
boot_Mdl.Const      = boot_Const;
boot_Mdl.A          = boot_A;
boot_Mdl.Sigma_eta  = boot_Sigma_eta;
boot_Mdl.comp_matr  = boot_Comp_Mat;
boot_Mdl.error  = boot_error;
boot_Mdl.y      = boot_y;

D = dupl(n);
%% MBB estimate of V
meanSigma_eta = vech(mean(boot_Sigma_eta,3)); 

SigSigmaEta = ... compute outer products for all bootstrap iterations
    arrayfun(@(b) (vech(boot_Sigma_eta(:,:,b))-meanSigma_eta)*(vech(boot_Sigma_eta(:,:,b))-meanSigma_eta)', 1:n_bootstrap,'UniformOutput', false);
V_MBB = ... MBB estimator of the covariance of Sigma_eta
    mean(reshape(cell2mat(SigSigmaEta),n*(n+1)/2,n*(n+1)/2,n_bootstrap),3);

V_MBB = T*V_MBB;
%% Bootstrapped estimate of V 
%  for b=1:n_bootstrap
%     iSigEta = inv( boot_Sigma_eta(:,:,b));
%     SigEta = boot_Sigma_eta(:,:,b);
%     for t=1:T-p
%         Upsilon2(:,t)= D'*(vec(iSigEta) - vec(iSigEta*boot_error(t,:,b)'*boot_error(t,:,b)*iSigEta));
%     end
%     iJ22 = pinv(full(D))*kron(SigEta,SigEta)*pinv(full(D))';
%     B = zeros(size(iJ22,1));
%     for h=0:T-p-1
%         B = B + ((h~=0)+1)*Upsilon2(:,1:end-h) * Upsilon2(:,h+1:end)';
%     end
%     V(:,:,b) = iJ22*B*iJ22;
% end
% Vhat = mean(V,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




function [b_Const, b_A, b_Sigma_eta, b_Comp_Mat, b_error, y_star] = MBB_iteration(u_blocks,u_center,T,n,p,y,Const,A,Mdl,block_length,numblocks,ratio,index)
u_tmp = zeros((block_length)*ratio,n);
for i=1:ratio
%     index = randi(numblocks);
    u_tmp((i-1)*block_length+1:i*block_length,:) = u_blocks(:,:,index(:,i));
end
u_tmp  = u_tmp - u_center; % recenter
%
u_star = u_tmp(1:end-(ratio*block_length-T),:);     % truncate

y_star = u_star;
y_star(1:p,:) = y(1:p,:); % initialize

for t=p+1:T % simulate
    y_star(t,:) = y_star(t,:) + Const';
    for l=1:p
        y_star(t,:) = y_star(t,:)' + ...
            A{l}*y_star(t-l,:)';
    end
end

[b_Mdl, se, ll, b_Mdl_errors] = estimate(Mdl,y_star); % estimation

b_Const = b_Mdl.Constant;
AL=[];
for l=1:p
    b_A(:,:,l) = b_Mdl.AR{1,l};
    AL = [AL b_A(:,:,l)];
end
b_Sigma_eta = b_Mdl.Covariance;
b_Comp_Mat = [AL;  eye(n*(p-1)) zeros(n*(p-1),n)];
b_error = b_Mdl_errors;
end