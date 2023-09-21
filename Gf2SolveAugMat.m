function [M, rank, EmptyRow, AugSol] = Gf2SolveAugMat(A)
[Nr Nc]=size(A);
Choosed=zeros(1,Nr); 
rank=0;

for j=1:Nc
  NewRankFlag=0 ;
  for index=1:Nr
    if( A(index,j)==1 & Choosed(index)==0) % find the row "index" of the first 1 from the top in column j 
      Choosed(index)=1;
      NewRankFlag=1;
      rank=rank+1;
      break;
    end
  end
  if NewRankFlag==1
    for Ir=1:index-1
      if( A(Ir,j) == 1)
        A(Ir,:) = xor( A(Ir,:), A(index,:)) ;
      end
    end 
    for Ir=index+1:Nr
      if( A(Ir,j) == 1)
        A(Ir,:) = xor( A(Ir,:), A(index,:)) ;
      end
    end
  end
end
  M=A;
  EmptyRow=find(Choosed-1);

AugSol=zeros(1,Nr) ;
for i=1:Nr
  if(Choosed(i)==1)
    AugSol(min(find(M(i,:)==1))) = M(i,end) ;
  end
end
  
  