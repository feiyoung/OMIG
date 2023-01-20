function [hX, hD, hHm, hBm, history] = OrMIG1(Xmis, group, type, q, dropout, eps2, maxIter, output)
% The function to conduct MIG by seperately exerting the
% identifiability conditions for H and (mu, B).
% Created by Wei Liu on 2020/01/18.
% Copyright @ 2019 Wei Liu. All rights reserved.
% dropout=0;
if(~exist('output', 'var'))
    output = 0;
end
ind_set = unique(group);
ng = length(ind_set);
if ~isempty(setdiff(1:ng, ind_set))
    error('Id number of types must match type!')
end
if ng ~= size(type,1)
    error('The number of groups must match with length of type!');
end

gcell = cell(1, ng);

for j = 1:ng
    gcell{j} = find(group==j);
end
[n,p] = size(Xmis);
omega = p^(-1);
% initialize
warning('off');
Xm = Xmis; Xm(isnan(Xm))=0;
hH = factorm(Xm, q, 0);
clear Xm;
eps1 = 1e-4;
hB = 0;
dBm = Inf; dH = Inf; dOmega = Inf;dc =Inf;
dOmega = max([dBm, dH]);
tmpB = zeros(p,q+1);tmpH = hH; tmpc = 1e7;
% maxIter = 50;
k = 1;
tic;
while k <= maxIter && dOmega > eps1 && dc > eps2
    hhB = [];
    for j = 1:ng
        
        [B1] = localupdateB2(Xmis, gcell{j}, hH, type(j,:));
        hhB = [hhB, B1];
    end
    hmu = hhB(1,:)';
    hB = hhB(2:end,:)';
    % ensure indentifiability.
    [B0tmp, ~] = qr(hB, 0);
    B0= B0tmp * diag(sort(sqrt(eig(hB'*hB)), 'descend'));
    
    sB = sign(B0(1,:));
    hB = B0.*repmat(sB,p,1); % ensure B first nonzero is positive
    %hB(1:4,:), B(1:4,:)
    hBm = [hmu, hB];
    dB = norm(hBm - tmpB, 'fro')/norm(hBm, 'fro');
    tmpB = hBm;
    if output
        fprintf('-------------------------------------------\n')
        fprintf('---------- B updation is finished!---------\n')
    end
    % given B^(1), update H^(1)
    H4 = localupdateH2(Xmis, gcell, hBm, type, dropout);
    hH0 = H4(:,2:end);
    
    [H0, ~] = qr(hH0, 0);
    hH1 = H0 * sqrt(n);
    sH0 = sign(hH0(1,:)).* sign(hH1(1,:));
    hH = hH1.* repmat(sH0,n,1);
    dH = norm(hH-tmpH, 'fro')/norm(hH, 'fro');
    tmpH = hH;
    if output
        fprintf('-------------------------------------------\n')
        fprintf('---------- H updation is finished!---------\n')
        fprintf('-------------------------------------------\n')
    end 
    hHm = [ones(n,1), hH];
    dOmega = max([dB, dH]);
    c = objMisfun(hHm, hBm, Xmis, omega, gcell, type);
    dc = abs(c - tmpc)/abs(tmpc);
    tmpc = c;
    if output
        fprintf('Iter %d \n', k);
        fprintf('dB= %4f, dH= %4f,dc= %4f, c=%4f \n', dB, dH,dc, c);
    end
    history.dB(k) = dB; history.dH(k) = dH; history.dc(k)=dc; history.c(k)=c;history.realIter = k;
    k = k+1;
end
history.maxIter = maxIter;
hD = hHm * hBm';
hX = Xmis;
for j= 1:ng
    switch type{j,1}
    case 'normal'
       hX(:, gcell{j}) = hD(:, gcell{j});
    case 'poisson'
       hX(:, gcell{j}) = round(exp(hD(:, gcell{j})));
    case 'binomial'
       hX(:, gcell{j}) = (1 ./ (exp(-hD(:, gcell{j})) +1) > 0.5) + 0;
    end
end
O = ~isnan(Xmis);
hX(O) = Xmis(O);