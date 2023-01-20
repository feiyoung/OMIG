 function  Uj = getScore(Xmisj, Hm, Bmj, typej, rb_vecj)
  % Xmisj=Xmisb(:,j); Hm=Hmb; Bmj= Bms(j,:); typej=type{group(j),1}; rb_vecj =rb_vec(j);
     n = length(Xmisj);
     o = (~isnan(Xmisj));
     Xmisj(~o) = 0;
     if (isempty(rb_vecj))
         OdQ = o/ mean(o);
     else
        OdQ = o / rb_vecj; 
     end
     
     
     switch typej
         case 'normal'
             
             muj = Hm * Bmj';
         case 'poisson'
             muj = exp(Hm * Bmj');
         case 'binomial'
             muj = 1 / (1+exp(-Hm * Bmj'));
     end
     % size(Xmisj), size(muj)
     Uj = Hm' * ((Xmisj - muj) .* OdQ);
     