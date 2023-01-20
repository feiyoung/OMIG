function [mea] = NAE(hX, X, Xmis, group)

hmu = mean(Xmis, 'omitnan');
Mu = repmat(hmu, size(Xmis,1), 1);
d = length(unique(group));
mea = zeros(1, d);
for i = 1:d
   mind1 = isnan(Xmis(:, group==i));
   hXi = hX(:, group==i);
   Xi = X(:, group==i);
   Mui = Mu(:, group==i);
   mea(i) = sum(abs(hXi(mind1) - Xi(mind1))) / sum(abs(Mui(mind1) - Xi(mind1)));
end

end