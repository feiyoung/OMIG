function [hX] = imputeFun(hD, Xmis, type, group)

ind_set = unique(group);
ng = length(ind_set);
gcell = cell(1, ng);

for j = 1:ng
  gcell{j} = find(group==j); 
end
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