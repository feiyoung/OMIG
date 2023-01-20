function [hesList] = getHessian(Xmis1, Hm, Bm, type, group, rb_vec)
% Xmis1 = Xmisb; Bm = Bms;
[n, p] = size(Xmis1);
q1 = size(Hm, 2);
O = (~isnan(Xmis1));
Xmis1(~O) = 0;
ng = size(type,1);
if(isempty(rb_vec))
    OdQ = O ./ repmat(mean(O)', n, 1);
else
   OdQ =  O./ repmat(rb_vec, n, 1);
end

gcell = cell(1,ng);
for j = 1:ng
   g1 = find(group ==j);
   gcell{j} = g1;
end

muMat = zeros(n, p);
for j = 1:ng
    if strcmp(type{j,1}, 'normal')
      muMat(:,gcell{j}) = Hm * Bm(gcell{j},:)';
    elseif strcmp(type{j,1}, 'poisson')
      muMat(:,gcell{j}) = exp(Hm * Bm(gcell{j},:)');  
    elseif strcmp(type{j,1}, 'binomial')
      muMat(:,gcell{j}) =1 / (1+exp(-Hm * Bm(gcell{j},:)'));
    else
       error('Unsupported variable type!');
    end
end
hesList = zeros(q1, q1, p);
for j = 1:p
    % j = 1;
    if strcmp(type{group(j),1}, 'normal')
     hesList(:,:,j) = Hm' * diag(OdQ(:,j)) * Hm;
    elseif strcmp(type{group(j),1}, 'poisson')
      hesList(:,:,j) = Hm' * diag(muMat(:,j) .* OdQ(:,j)) * Hm;
    elseif strcmp(type{group(j),1}, 'binomial')
      hesList(:,:,j) = Hm' * diag(muMat(:,j) .* (1-muMat(:,j)) .* OdQ(:,j)) * Hm;
    else
       error('Unsupported variable type!');
    end
end