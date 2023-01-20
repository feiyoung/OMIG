function idxList = batchGroup(n, nbatch1, nbatch, nleast)

n1 = nbatch;
if nbatch1 ==0
   J = floor(n/n1);
   idxList = cell(1,J);
   for j = 1:J
       idxList{j} = ((j-1)*n1+1):(j*n1);
   end
   if (n-n1*J >= nleast)
      J = J+1;
      idxList{J} = ((j*n1)+1):n;
   else
       idxList{J} = [idxList{J}, ((j*n1)+1):n];
   end
elseif (nbatch1>0)
   
   n2 = n - nbatch1;
   J = floor(n2/n1)+1;
   idxList = cell(1,J);
   idxList{1} = 1:nbatch1;
   for j = 2:J
      idxList{j} = nbatch1 +  (((j-2)*n1+1):((j-1)*n1));
   end
   if((n2 -n1*J) >= nleast)
       J = J+1;
       idxList{J} = ((j*n1)+1):n;
   elseif(nbatch1+ ((j-1)*n1)<n)
      idxList{J} = [idxList{J}, (nbatch1 + ((j-1)*n1+1):n)];
   end
end
end