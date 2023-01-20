function Hm = updateHs(Xmis, Bm, Hm, gcell, type, rb_vec, parallel, lambda)
[n, p] = size(Xmis);
O = ~isnan(Xmis); % index of observation
if(isempty(rb_vec))
   OdQ = O ./ repmat(mean(O), n, 1); % o_{ij} / q_j.
else
   OdQ = O ./ repmat(rb_vec, n, 1); 
end
Xmis(isnan(Xmis))=0;
B = Bm(:,2:end);
q = size(Hm, 2) - 1;
ng = size(type,1);

% mean matrix
mucell = cell(1,ng);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        mucell{j} = Hm * Bm(gcell{j},:)';
    case 'poisson'
        mucell{j} = exp(Hm * Bm(gcell{j},:)');
    case 'binomial'
        jvec = gcell{j};
        Xj = Xmis(O(:,jvec(1))==1, jvec(1));
        ntrail_j = length(unique(Xj))-1;
        mucell{j} =ntrail_j*1./(1 + exp(-Hm * Bm(gcell{j},:)'));
    end
end
% mucell{3}(1:4, 1:8)
% socre matrix
df2 = zeros(n, q);
for j = 1:ng
    df2 = df2 + ((Xmis(:, gcell{j})- mucell{j}).* OdQ(:, gcell{j}) )* B(gcell{j},:);
end
df2  =  df2 - lambda* Hm(:, 2:end);
% % Hessian matrix or information matrix
d2f = cell(n,1);
if parallel
    parfor i = 1:n
        d2f{i} = hessianMat(i, B, gcell,mucell, OdQ(i,:), type, lambda);
    end
else
    for i = 1:n
        d2f{i} = hessianMat(i, B, gcell,mucell, OdQ(i,:), type, lambda);
    end
end


cell1 = cell(n,1);
if parallel
    parfor i = 1:n
        cell1{i} = df2(i,:)';
    end
else
    for i = 1:n
    cell1{i} = df2(i,:)';
end
end

H2update = cellfun(@(x, y) (x+ 1e-20*eye(q))\y, d2f, cell1, 'UniformOutput', 0);
% H2 = Hm - cell2array(H2update)';
H2updatemat = reshape(cell2mat(H2update), q, n);
Hm(:,2:end) = Hm(:, 2:end) + H2updatemat';


function [d2f] = hessianMat(i, B, gcell,mucell, OdQi, type, lambda)
   q = size(B,2);
   ng = size(type,1);
   Bng = zeros(q,q);
   for j = 1:ng
        switch type{j,1}
    case 'normal'
        %W = diag(1./ (std(X(:,gcell{j})).^2));
        Bng = Bng + B(gcell{j},:)'* diag(OdQi(gcell{j}))*B(gcell{j},:);
    case 'poisson'
        Bng = Bng + B(gcell{j},:)'* diag(mucell{j}(i,:).* OdQi(gcell{j})) * B(gcell{j},:);
        %Bng = Bng + (repmat(mucell{j}(i,:), q, 1)'.* B(gcell{j},:))'* B(gcell{j},:);
    case 'binomial'
        %Bng = Bng + (repmat(mucell{j}(i,:), q,1)' .* B(gcell{j},:))' *(repmat(1-mucell{j}(i,:), q,1)' .* B(gcell{j},:));
        Bng = Bng + (B(gcell{j},:))' * diag(mucell{j}(i,:).*(1-mucell{j}(i,:)).* OdQi(gcell{j})) * B(gcell{j},:);

        end
   end
    d2f = Bng + lambda* eye(q,q);