function [X_imp, hXall, hD, Hms, Bms, time_vec, NAE_mat, AE_mat, corB_vec, corH_vec, errD] = OMIG(Xmis, group, type, q, idxList, lambda,verbose,parallel, X, Bm0, H0)


J = length(idxList);
time_vec = zeros(1,J);
errD = [];
if(isempty(X))
   NAE_mat = []; 
   AE_mat = [];
else
   NAE_mat = zeros(J, size(type,1)); 
   AE_mat = zeros(J, size(type,1)); 
end
if(isempty(Bm0))
    corB_vec = [];
else
   corB_vec = zeros(1,J);
   errD = [];
end

if(isempty(H0))
    corH_vec = [];
else
   corH_vec = zeros(1,J);
   errD = [];
end
X_imp = Xmis;

b = 1;
[n, p] = size(Xmis);
if((~isempty(H0)) && (~isempty(Bm0)))
   errD = zeros(1,J);
   D0 = [ones(n,1), H0] * Bm0';
end

Xmisb = Xmis(idxList{b},:);
ns = size(Xmisb, 1);
nrb_vec = size(Xmisb, 1) * mean(~isnan(Xmisb));
rb_vec = nrb_vec / ns;
tic; 
try
  [hXb, hD, Hm, Bms] = OrMIG(Xmisb, group, type, q, 1e-4,...
    20,lambda, verbose, parallel); %% here require to change 5 to 20
  X_imp(idxList{b},:) = hXb;
catch
  hmu = mean(Xmis, 'omitnan'); 
  X0 = Xmisb - repmat(hmu', ns, 1);
  X0(isnan(Xmisb)) = 0;
  [hatH, hatB]  = factorm(X0, q);
  Hm = [ones(ns, 1) hatH]; Bms = [hmu', hatB];
  hXb = Hm * Bms';
end
time1 = toc;
time_vec(b) = time1;
if(~isempty(X))
   NAE_mat(b,:) = NAE(hXb, X(idxList{b},:), Xmisb, group); 
   AE_mat(b,:) = AE(hXb, X(idxList{b},:), Xmisb, group); 
end
if(~isempty(Bm0))
  [~,~,rtmp] = canoncorr(Bms, Bm0);
  corB_vec(b) = rtmp(q+1);
end
if(~isempty(H0))
  [~,~,rtmp] = canoncorr(Hm(:,2:end), H0(idxList{b},:));
  corH_vec(b) = rtmp(q);
end

Hms = Hm;
hesLists = getHessian(Xmisb, Hms, Bms, type, group, rb_vec);
hXb = imputeFun(Hms * Bms', Xmisb, type, group);
X_imp(idxList{b},:) = hXb;

if((~isempty(H0)) && (~isempty(Bm0)))
   idx = 1:idxList{b}(end);
   errD(b) = sum(sum(abs(Hms*Bms' - D0(idx,:)))) / numel(D0(idx,:));
end

for b = 2:J
    % b = 2;
    idx = idxList{b};
    Xmisb = Xmis(idx,:);
    nrb_vec = nrb_vec + size(Xmisb,1) * mean(~isnan(Xmisb));
    ns = size(Xmisb, 1) + ns;
    rb_vec = nrb_vec / ns;
    
    ind_set = unique(group);
    ng = length(ind_set);
    gcell = cell(1, ng);
    for j = 1:ng
       g1 = find(group ==j);
       gcell{j} = g1;
       if strcmp(type{j,1}, 'binomial')
          N = max(Xmisb(:, g1),  'omitnan');
          Xmisb(:, g1) = Xmisb(:, g1) / N;
       end
    end
    
    tic;
    hmu = mean(Xmisb, 'omitnan');
    X0 = Xmisb - repmat(hmu, size(Xmisb,1), 1);
    X0(isnan(Xmisb)) = 0;
    [hatH, ~]  = factorm(X0, q);
    Hmb = [ones(size(X0,1), 1) hatH]; 
    maxIter = 15;
    negloglike_seq = zeros(maxIter,1);
    negloglike_seq(1) = 1e10;
    for iter = 2:maxIter
        % iter = 2;
       try
           Hmb = updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec, parallel, lambda);
       catch
           Hmb = Hmb; 
       end
       negloglike_seq(iter) = objMisfun(Hmb, Bms, Xmisb, gcell, type);
       dc = abs(negloglike_seq(iter) - negloglike_seq(iter-1))/abs(negloglike_seq(iter-1));
       fprintf('dc = %4f \n', dc);
       if(dc < 1e-9) 
           break;
       end
    end
    
    %Given Hm, update Upsilon
    hesList2 = getHessian(Xmisb, Hmb, Bms, type, group, rb_vec);
    hesLists1 = hesLists;
    for j = 1:p
        hesLists(:,:,j) = hesLists1(:,:,j) + hesList2(:,:,j);
    end
    maxIter = 15;
    negloglike_seq = zeros(maxIter,1);
    negloglike_seq(1) = 1e10;
    Bm1 = Bms;
    for iter = 2 : maxIter
        % iter = 2;
        for j = 1:p
            %j = 1;
            Ubj = getScore(Xmisb(:,j), Hmb, Bms(j,:), type{group(j),1}, rb_vec(j));
            Uj = hesLists1(:,:,j) * (Bm1(j,:) - Bms(j,:))' + Ubj;
            
            Bms(j,:) = Bms(j,:) + Uj'/(hesLists(:,:,j));
        end
        n_nan = sum(sum(isnan(Bms)));
        Bms(isnan(Bms)) = rand(n_nan,1);
        negloglike_seq(iter) = objMisfun(Hmb, Bms, Xmisb, gcell, type);
       dc = abs(negloglike_seq(iter) - negloglike_seq(iter-1))/abs(negloglike_seq(iter-1));
       fprintf('dc = %4f \n', dc);
       if(dc < 1e-9) 
           break;
       end
        
    end
    
    % re-update H
    maxIter = 15;
    negloglike_seq = zeros(maxIter,1);
    negloglike_seq(1) = 1e10;
    for iter = 2:maxIter
        % iter = 2;
       try
           Hmb = updateHs(Xmisb, Bms, Hmb, gcell, type, rb_vec, parallel, lambda);
       catch
           Hmb = Hmb;
       end
       negloglike_seq(iter) = objMisfun(Hmb, Bms, Xmisb, gcell, type);
       dc = abs(negloglike_seq(iter) - negloglike_seq(iter-1))/abs(negloglike_seq(iter-1));
       fprintf('dc = %4f \n', dc);
       if(dc < 1e-9) 
           break;
       end
    end
    Hms = [Hms; Hmb];
    hXb = imputeFun(Hmb * Bms', Xmisb, type, group);
    X_imp(idxList{b},:) = hXb;
     time1 = toc;
     time_vec(b) = time1;
    if(~isempty(X))
      NAE_mat(b,:) = NAE(hXb, X(idxList{b},:), Xmisb, group); 
      AE_mat(b,:) =  AE(hXb, X(idxList{b},:), Xmisb, group); 
    end
    if(~isempty(Bm0))
      [~,~,rtmp] = canoncorr(Bms, Bm0);
      corB_vec(b) = rtmp(q+1);
    end
    if(~isempty(H0))
      [~,~,rtmp] = canoncorr(Hmb(:,2:end), H0(idxList{b},:));
      corH_vec(b) = rtmp(q);
    end
    if((~isempty(H0)) && (~isempty(Bm0)))
       idx = 1:idxList{b}(end);
       errD(b) = sum(sum(abs(Hms*Bms' - D0(idx,:)))) / numel(D0(idx,:));
    end
    
  
end
hXall = imputeFun(Hms * Bms', Xmis, type, group);
    
