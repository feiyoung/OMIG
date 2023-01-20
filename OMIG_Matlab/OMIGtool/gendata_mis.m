function [Xmis, X, H, Bm, group] = gendata_mis(seed, n, p, vtype, q, mis_vec, rho, mechanism)
% type = {'homonorm', 'heternorm', 'pois', 'norm_pois', 'pois_bino', 'npb'}
if(~exist('n', 'var'))
    n = 300;
end
if(~exist('p', 'var'))
    p = 50;
end
if(~exist('vtype', 'var') || isempty(vtype))
     vtype = 'homonorm';
end
if(~exist('q', 'var'))
    q = 6;
end
if(~exist('mis_vec', 'var'))
    mis_vec = 0.3;
end
if(~exist('rho', 'var'))
    rho = 3;
end

if(~exist('mechanism', 'var'))
    mechanism = 'MCAR';
end

if(length(rho)==1) 
    rho = [rho, 1.5];
end
factor_term = rho(2);
rng(1); % to fix B and mu
        Z =  randn(p,q);
        [B1, ~,~] = qr(Z, 0);
        B = rho(1)* B1 * diag(sqrt(sort(eig(Z'*Z),'descend'))) / sqrt(n); % sort the eigenvectors by decreasing eigen values.
        % V = B'*B/p;
        mu = 0.4*randn(1,p);
        rng(seed);  % For reproducibility
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        
switch vtype
    case 'homonorm'
     
        X = repmat(mu, n, 1) + H * B' + mvnrnd(zeros(1,p), diag(ones(p,1)), n);
        group = ones(1, p);
    case 'heternorm'
        
        sigmas = 0.1+2*rand(p,1); % heteroskedasticity
        X = repmat(mu, n, 1) + H * B' + mvnrnd(zeros(1,p), diag(sigmas), n);
        group = ones(1, p);
    case 'pois'
       
        % genarate X
        g2 = 1:p; % Poisson: exp
        %  /max(B0[g1,])* factor_term
        mu2 = exp(H* B(g2,:)'+ repmat(mu(g2), n, 1)); % poisson distribution
        X = poissrnd(mu2);
        group = ones(1, p);
    case 'norm_pois'
        % genarate X
        g1 = 1:floor(p/2); % identity, normal
        g2 = (floor(p/2)+1):p; % Poisson: exp
        
        mu1 = H* B(g1,:)' + repmat(mu(g1), n, 1); % normal dstribution.
        X1 = normrnd(mu1, 0.01);
        B(g2,:) = B(g2,:) / max(max(B(g2,:))) * factor_term; 
        mu2 = exp(H* B(g2,:)'+ repmat(mu(g2), n, 1)); % poisson distribution
        X2 = poissrnd(mu2);
        
        X = [X1 X2];
        group = [ones(1, length(g1)) 2*ones(1, length(g2))];
    case 'pois_bino'
        % genarate X
        g1 = 1:floor(p/2); % poisson
        g2 = (floor(p/2)+1):p; % Poisson: exp
        
        B(g1,:) = B(g1,:) / max(max(B(g1,:))) * factor_term;
        mu1 = exp(H* B(g1,:)'+repmat(mu(g1), n, 1)); % poisson distribution
        X1 = poissrnd(mu1);
        
        mu2 = 1./(1 + exp(-H*B(g2,:)'- repmat(mu(g2), n, 1))); % binary distribution.
        
        X2 = binornd(1, mu2);
        X = [X1 X2];
        group = [ones(1, length(g1)) 2*ones(1, length(g2))];
    case 'npb'
        % genarate X
        g1 = 1:floor(p/3); % identity, normal
        g2 = (floor(p/3)+1):floor(2*p/3); % Poisson: exp
        g3 = (floor(2*p/3) + 1):p; % Bernoulli
        
        mu1 = H* B(g1,:)' + repmat(mu(g1), n, 1); % normal dstribution.
        X1 = normrnd(mu1, 1);
        B(g2,:) = B(g2,:) / max(max(B(g2,:))) * factor_term;
        mu2 = exp(H* B(g2,:)' + repmat(mu(g2), n, 1)); % poisson distribution
        X2 = poissrnd(mu2);
        
        mu3 = 1./(1 + exp(-H*B(g3,:)' - repmat(mu(g3), n, 1))); % binary distribution.
        X3 = binornd(1, mu3);
        
        X = [X1 X2 X3];
        group = [ones(1, length(g1)) 2*ones(1, length(g2)), 3*ones(1, length(g3))];
    case  'bino'
        g = 1:p;
        mu1 = 1./(1 + exp(-H*B(g,:)'- repmat(mu, n, 1))); % binary distribution.
        N = 1;
        X = binornd(N, mu1);
        group = ones(1, p);
    case 'norm_bino'
        g1 = 1:floor(p/2); % normal
        g2 = (floor(p/2)+1):p; %  binomial
        
        mu1 = H* B(g1,:)'+ repmat(mu(g1), n, 1); % normal dstribution.
        X1 = normrnd(mu1, 1);
        mu2 = 1./(1 + exp(-H*B(g2,:)' - repmat(mu(g2), n, 1))); % binary distribution.
        X2 = binornd(1, mu2);
        X = [X1 X2];
         group = [ones(1, length(g1)) 2*ones(1, length(g2))];
end
Bm = [mu', B];
% generate missing data
Xmis = X;
hasBino = 0;
if any(strcmp({'npb', 'pois_bino', 'bino'},vtype))
    hasBino = 1;
    switch vtype
        case 'bino'
           idx_bino  = 1:p;
        case 'pois_bino'
            idx_bino = (floor(p/2)+1):p;
        case 'npb'
            idx_bino = (floor(2*p/3)+1):p;
    end
end
if strcmp(mechanism, 'MCAR')
    n_mis = length(mis_vec);
    for kk = 1:n_mis
        for j = 1:p
            if rem(j, n_mis) == (kk-1)
               mis_rate = mis_vec(kk);
               misj_ind = datasample(1:n, floor(n*mis_rate));
               Xmis(misj_ind, j) = NaN;
            end
        end
    end
elseif strcmp(mechanism, 'MNAR')
         n_mis = length(mis_vec);
    for kk = 1:n_mis
        for j = 1:p
            if rem(j, n_mis) == (kk-1)
               mis_rate = mis_vec(kk);
               if hasBino && isempty(setdiff(j,idx_bino))
                   misj_ind = datasample(1:n, floor(n*mis_rate));
               else
                   [~,misj_ind] = sort(X(:,j), 'descend');
               end
               
               Xmis(misj_ind(1:floor(n*mis_rate)), j) = NaN;
            end
        end
    end
elseif strcmp(mechanism, 'MAR')  
      n_mis = length(mis_vec);
      p1 = floor(p/2);
      XX = [X(:, (p1+1):p), X(:, 1:p1)];
    for kk = 1:n_mis
        for j = 1:p
            if rem(j, n_mis) == (kk-1)
               mis_rate = mis_vec(kk);
               if j <= floor(p/2)
                   [~,misj_ind] = sort(XX(:,j), 'descend');
                   misj_ind = misj_ind(1:floor(n*mis_rate));
                   Xmis(misj_ind, j) = NaN;
               else
                   
                   misj_ind = datasample(1:n, floor(n*mis_rate));
                   Xmis(misj_ind, j) = NaN;
               end
            end
        end
    end 
end
end