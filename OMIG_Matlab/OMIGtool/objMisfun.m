function obj = objMisfun(Hm, Bm, Xmis,  gcell, type)
% function to evaluate the value of the objective function (the conditional
% loglikelihood function).
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         June. 25, 2021
% Copyright (c) 2021, Liu Wei
% All rights reserved.

n = size(Xmis, 1);
O = ~isnan(Xmis); % index of observation
OdQ = O ; % ./ repmat(mean(O), n, 1); % o_{ij} / q_j.
Xmis(isnan(Xmis))=0;
eps1 = 1e-20;
Bh = Hm * Bm';
ng = size(type,1);
Q = zeros(size(Xmis));
for j = 1:ng
    switch type{j,1}
    case 'normal'        
        Q(:, gcell{j}) = ((Xmis(:, gcell{j})- Bh(:, gcell{j})) .* OdQ(:, gcell{j}) ).^2;
    case 'poisson'
        me = exp(Bh(:,gcell{j}));
        Q(:, gcell{j}) = -(log(poisspdf(Xmis(:, gcell{j}), me)+eps1) .* OdQ(:, gcell{j}));
    case 'binomial'
        me3 = 1 ./(1+exp(-Bh(:,gcell{j})));
        Q(:,gcell{j}) = -(Xmis(:,gcell{j}).*log(me3+eps1) + (1-Xmis(:,gcell{j})).* log(1-me3+eps1)) .* OdQ(:, gcell{j}) ;
    end
end

obj = mean(Q(O));
