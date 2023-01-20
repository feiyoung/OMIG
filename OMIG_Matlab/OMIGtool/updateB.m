function Bm2 = updateB(Xmis, hBm, hHm, gcell, type)
% function to conduct one-step updating for efficiently estimate latent factor matrix H.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.
O = ~isnan(Xmis); % index of observation
Xmis(isnan(Xmis))=0;
[~, p] = size(Xmis);
q = size(hHm, 2) - 1;
ng = size(type,1);

% mean matrix
Bm2 = zeros(p, q+1);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        mutypej = hHm * hBm(gcell{j},:)'; % n*p1
        scorej = hHm' * ((Xmis(:, gcell{j}) -mutypej).* O(:, gcell{j})); % (q+1) * p1
        Hesj =  hHm'* hHm;
        Bm2(gcell{j},:) = hBm(gcell{j},:) + scorej' / Hesj;
    case 'poisson'
        jvec = gcell{j};
        mutypej = exp(hHm * hBm(jvec,:)');
        scorej = hHm' * ((Xmis(:, jvec) -mutypej).* O(:, gcell{j})); % (q+1) * p1
        p1 = length(jvec);
        for jl = 1:p1
            Hestypejjl =  hHm'* diag(mutypej(:,jl).* O(:, jvec(jl))) * hHm; % take the observed index to Hessian.
            Bm2(jvec(jl),:) = hBm(jvec(jl),:) + scorej(:,j)' / Hestypejjl;
        end
        
    case 'binomial'
        jvec = gcell{j};
        Xj = Xmis(:, jvec);
        ntrail_j = length(unique(Xj(:,1)))-1;
        mutypej =ntrail_j*1./(1 + exp(-hHm * hBm(jvec,:)'));
        scorej = hHm' * ((Xj -mutypej).*  O(:, gcell{j})); % (q+1) * p1
        p1 = length(jvec);
        
        for jl = 1:p1
            Hestypejjl =  hHm'* diag(mutypej(:,jl).*(1-mutypej(:,jl)))  * hHm;
            Bm2(jvec(jl),:) = hBm(jvec(jl),:) + scorej(:,j)' / Hestypejjl;
        end
    end
end

