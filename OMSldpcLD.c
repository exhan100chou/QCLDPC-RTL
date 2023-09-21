#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CirculantSize 256
#define Weight 4
#define InitialMatrixSeed 123
#define LayerSize 256
#define RateParameter 10 // 8 code rate = 1 - (1/RateParameter) ;

#define BitNodeNum (CirculantSize*RateParameter*Weight)  
#define CheckNodeNum (CirculantSize*Weight)
#define ShortenNum 0

#define OffsetFactor 0.75

#define BitNodeDeg Weight
#define CheckNodeDeg (Weight*RateParameter)
#define MaxIterNum 8

#define Inf 1000000
int             recd[1048576], data[1048576], bb[548576];
int             seed;
int GLeadingvector[CheckNodeNum*(BitNodeNum-CheckNodeNum)];
char H[CheckNodeNum*BitNodeNum];
char D[CheckNodeNum*CheckNodeNum];
char AugD[CheckNodeNum*(CheckNodeNum+1)];
char AugH[CheckNodeNum*(BitNodeNum+1)];
int di[CheckNodeNum+1]={0};
int diH[CheckNodeNum+1]={0};
int swaprow[10*CheckNodeNum]={0};int swaprowH[10*CheckNodeNum]={0};
int ZiT[CheckNodeNum];
typedef struct
{
  double ExtInf[BitNodeDeg];
  double RecLlr;
  unsigned char TempHD;
  int index;
}BitNode;

typedef struct
{
  int BitNodeIdx[CheckNodeDeg];
}CheckNode;

typedef struct
{
  BitNode BN[BitNodeNum];
  CheckNode CN[CheckNodeNum] ;
  int HqcLeader[CheckNodeNum/CirculantSize][BitNodeNum/CirculantSize];
  int GqcLeader[(BitNodeNum-CheckNodeNum)/CirculantSize][CheckNodeNum/CirculantSize];
}Ldpc;

int SignFunction(double x)
{
  return (x>=0) ? 1 : -1 ;
}

double PhiFunction(double x)
{
   double y = 0;
   y = (x>0) ? x : -x;
   if(y > 35)
     return 1.2e-015;
   else if(y < 1.2e-015)
     return 35;
   else
     return log((exp(y)+1)/(exp(y)-1));
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
/* Copy from Numerical Recipes in C  2nd edt  Page 279
  set idum to be any integer value
  Input a seed , this sub-function will return a random number between [0 ,1] */

{
  long k ;
  double ans;

  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum <0)  *idum+=IM;
  ans=(double)(AM*(*idum));
  return ans;
}

/*Gaussian_noise mean=0 , variance=1 */
double gaussian_noise(long *seed_ptr)
{
  double noise,s,v[2];
  do{
      v[0] = 2*ran0(seed_ptr) - 1;
	   v[1] = 2*ran0(seed_ptr) - 1;
    	s = v[0]*v[0] + v[1]*v[1];
    } while(s>=1);
  noise = v[0]*sqrt(-2*log(s)/s);
  return(noise);
}

void RandomGenHMatrix(Ldpc *QC436)
{
  int i,j;
  int k;
  long MatrixSeed = InitialMatrixSeed;
  FILE *out1;

  for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{
	  if(i!=0 && j!=0) // Not in the first column or the first row
	    QC436->HqcLeader[i][j] = (int) (ran0(&MatrixSeed) * CirculantSize) ;
 	  else
		QC436->HqcLeader[i][j] = 0 ;

	  for(k=0 ; k<CirculantSize ; k++)
	    QC436->CN[i*CirculantSize+k].BitNodeIdx[j] = j*CirculantSize+ ( (QC436->HqcLeader[i][j]+k)%CirculantSize ) ;
    }

  out1=fopen("RandomHMatrix.txt","a+t");
  for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
  {
    fprintf(out1, "\n"); 
    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{
       fprintf(out1, "%3d ", QC436->HqcLeader[i][j]); 
	}
  }
  fclose(out1);
}
int swapnum;
	int i,j,r,v,k,l,temp,w;int t, c;char tmp;
// Subfunction to swap row x and row y in Matrix 
void SwapRows(char *M, const int RowNum, const int ColNum, const int x, const int y) 
{
  for(i=0 ; i<ColNum ; i++) 
  {
    tmp = M[x*ColNum + i];
    M[x*ColNum + i] = M[y*ColNum + i];
    M[y*ColNum + i] = tmp;
  }// end of for
}

// Subfunction to substract all elements in row x by all elements in row y
void RawSubstract(char *M, const int RowNum, const int ColNum, const int x, const int y) 
{

  for(i = 0; i < ColNum; i++) 
    M[x*ColNum + i] ^= M[y*ColNum + i]; 
}

// RankComputationOfCirMatrix to compute the rank of a circulant matrix
// inputs are:
//   char *in  : input square matrix with single row representation
//   int n     : row or column size of matrix => Matrix size n*n
// output is the the number of dependent columns
//   int *di
void RankComputationOfCirMatrix(char *M, const int RowNum, const int ColNum, int *di,int *swaprow) 
{
  k=0;
    
  // ------Gaussian step: -----------//
  for(i=0 ; i<RowNum-1; i++) 
  {
    // first, we need to find a row that have "1" at location [i, i]
    if(M[i*ColNum + i] == 0) 
    {
      for(j = i+1; j < RowNum; j++) 
      {
        if(M[j*ColNum + i] == 1 ) 
        {
          // find the row with "1" at location [j, i]
          // then we swap row i and row j
          SwapRows(M, RowNum, ColNum, i, j);
		  swaprow[k]=i;swaprow[k+1]=j;
		  k=k+2;
          break;
        }// end of if 
      }// end of for 
      // if we can't find a row with "1" at location, 
      if(j == RowNum) 
	  {
       // printf("Column %d Fails\n", i);
		di[ i ] ++ ;
	  }// end of if
	}// end of if(in[i*(n+1)] == 0) 

    // and now we let all row j below the current row that have "1" at location [j, i]
    // to become 0 at location [j, i]
    for(j = i+1; j < RowNum; ++j) 
	{
      if(M[j*ColNum + i] == 1) 
	  {
        RawSubstract(M, RowNum, ColNum, j, i);
	  }//end of if(in[j*n+i] == 1) 
	}// end of for(j = i+1; j < n; ++j) 
  } // end of for (Gaussian step)

  // after doing n-1 lines, check the last element of the last row, whether it's "1"?
  if(M[RowNum*ColNum -1] == 0) 
  {
  //  printf("Column %d Fails\n", i);
	di[ i ] ++ ;
  }// end of if(in[n*n-1] == 0) 
  swapnum=k;
 
}
void RankComputationOfCirMatrixH(char *M, const int RowNum, const int ColNum, int *di,int *swaprow) 
{
  k=0;
    
  // ------Gaussian step: -----------//
  for(i=0 ; i<RowNum-1; i++) 
  {
    // first, we need to find a row that have "1" at location [i, i]
    if(M[i*ColNum + i + (ColNum-RowNum)] == 0) 
    {
      for(j = i+1; j < RowNum; j++) 
      {
        if(M[j*ColNum + i + (ColNum-RowNum)] == 1 ) 
        {
          // find the row with "1" at location [j, i]
          // then we swap row i and row j
          SwapRows(M, RowNum, ColNum, i, j);
		  swaprow[k]=i;swaprow[k+1]=j;
		  k=k+2;
          break;
        }// end of if 
      }// end of for 
      // if we can't find a row with "1" at location, 
      if(j == RowNum) 
	  {
       // printf("Column %d Fails\n", i);
		di[ i ] ++ ;
	  }// end of if
	}// end of if(in[i*(n+1)] == 0) 

    // and now we let all row j below the current row that have "1" at location [j, i]
    // to become 0 at location [j, i]
    for(j = i+1; j < RowNum; ++j) 
	{
      if(M[j*ColNum + i+ (ColNum-RowNum)] == 1) 
	  {
        RawSubstract(M, RowNum, ColNum, j, i);
	  }//end of if(in[j*n+i] == 1) 
	}// end of for(j = i+1; j < n; ++j) 
  } // end of for (Gaussian step)

  // after doing n-1 lines, check the last element of the last row, whether it's "1"?
  if(M[RowNum*ColNum -1] == 0) 
  {
  //  printf("Column %d Fails\n", i);
	di[ i ] ++ ;
  }// end of if(in[n*n-1] == 0) 
  swapnum=k;
 
}
void BackwardSubstitution(int *di,char *AugD,int *ZiT,char *D,int *swaprow)
{  
	
	
	//FILE *out1;

 // out1=fopen("di.txt","a+t");
  for(r=0 ; r<CheckNodeNum ; r++)
  {
	if(r%CirculantSize == 0)
	//  fprintf(out1, "\n"); 
    
//	fprintf(out1, "\n"); 
    for(v=0 ; v<(CheckNodeNum+1); v++)
	{
	  if(v%CirculantSize == 0)
	  //  fprintf(out1, " "); 
      //Hans modify: delete the failed column to simplify the backward substitution
        if(di[v]!=0 && AugD[r*(CheckNodeNum+1) + v]==1)AugD[r*(CheckNodeNum+1) + v]=0;
	    if(di[r]!=0 && AugD[r*(CheckNodeNum+1) + v]==1)AugD[r*(CheckNodeNum+1) + v]=0;
    //  fprintf(out1, "%d", AugD[r*(CheckNodeNum+1) + v]); 
	}
  }
 // fprintf(out1, "\n");
//  fclose(out1);


// To Hans: augmented matrix here
	//Hans modify:Implement backward substitution
	for(k=CheckNodeNum-1 ; k>=0 ; k--)
	if(AugD[k*(CheckNodeNum+1)+k]==1)
	{
		temp=0;
	  // use augmented matrix to solve D ZiT = MiuT
       if(AugD[ (k+1)*(CheckNodeNum+1)-1]==1)temp^=1;
	   for(l=k+1;l<CheckNodeNum;l++)
	   {
         
		 if(AugD[k*(CheckNodeNum+1)+l]==1 && ZiT[l]==1 && di[l]==0)
             temp^= 1;
	        ZiT[k]=(char)temp;
	   }
	}
	 for(r=swapnum-1 ; r>0; r=r-2)
	 { 
       SwapRows(AugD, CheckNodeNum, CheckNodeNum+1, swaprow[r-1], swaprow[r]);
	 }

	//Hans modify:Check the backward substitution work or not
/*	for(j=0 ; j<CheckNodeNum ; j++)
    {
       temp=0;
      if(di[j]==0 )
	  {
      for(k=0 ; k<CheckNodeNum ; k++)
	  {
         if(AugD[j*(CheckNodeNum+1)+k]==1 && ZiT[k]==1&& di[k]==0   )
           temp^=1;
	  }
	  if(temp!=AugD[(j)*(CheckNodeNum+1)+CheckNodeNum])
	  {
		  printf("backward substitution fail\n");
		  break;
	  }
      }
	}//End of for(j=0 ; j<CheckNodeNum ; j++)
	*/
}
void ComputeEncoder(Ldpc *QC440S,int *GLeadingvector)
{
  

 // FILE *out1;




  
  // GqcLeader[t-c][c] 
  t = BitNodeNum/CirculantSize ;
  c = CheckNodeNum/CirculantSize ;
   
  for(i=0 ; i<CheckNodeNum ; i++)
    for(j=0 ; j<BitNodeNum ; j++)
	{
	  H[i*BitNodeNum + j] = 0 ;
    }
  for(i=0 ; i<CheckNodeNum ; i++)
  {
	for(j=0 ; j<t-c ; j++)
      GLeadingvector[i+j*CheckNodeNum]= -1;
  }
  for(i=0 ; i<CheckNodeNum ; i++)
  {
	for(j=0 ; j<CheckNodeDeg ; j++)
	H[ i*BitNodeNum + QC440S->CN[i].BitNodeIdx[j] ] = 1 ;
  }
/*
  out1=fopen("H.txt","a+t");
  for(i=0 ; i<CheckNodeNum ; i++)
  {
	if(i%CirculantSize == 0)
	  fprintf(out1, "\n"); 
    
	fprintf(out1, "\n"); 
    for(j=0 ; j<BitNodeNum ; j++)
	{
	  if(j%CirculantSize == 0)
	    fprintf(out1, " "); 

      fprintf(out1, "%d", H[i*BitNodeNum + j]); 
	}
  }
  fclose(out1);*/

// To Hans: should modified the subfunction to compute Rank of H and store in diH
//          I am not sure if row exchange should NOT be used?
  //RankComputationOfCirMatrix(H, CheckNodeNum, BitNodeNum, diH);

 /* out1=fopen("diH.txt","w");
  for(i=0 ; i<Weight*CirculantSize ; i++)
  {
    fprintf(out1, "%d ", diH[i]); 
  }
  fclose(out1);*/

  for(i=0 ; i<CheckNodeNum ; i++)
    for(j=0 ; j<CheckNodeNum ; j++)
	{
	  D[i*CheckNodeNum+j] = 0 ;
    }

// To Hans: Here I did not make sure that the last (rightest) part of H has the same rank with H
//          if not, we need to exchange circulant of column 
  for(i=0 ; i<CheckNodeNum ; i++)
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  if(QC440S->CN[i].BitNodeIdx[j] >= BitNodeNum-CheckNodeNum)
		D[i*CheckNodeNum + QC440S->CN[i].BitNodeIdx[j]-(BitNodeNum-CheckNodeNum)] = 1 ;
	}
  
 // out1=fopen("D.txt","a+t");
/*  for(i=0 ; i<CheckNodeNum ; i++)
  {
	if(i%CirculantSize == 0)
	  fprintf(out1, "\n"); 
    
	fprintf(out1, "\n"); 
    for(j=0 ; j<CheckNodeNum ; j++)
	{
	  if(j%CirculantSize == 0)
	    fprintf(out1, " "); 

      fprintf(out1, "%d", D[i*CheckNodeNum + j]); 
	}
  }
  fclose(out1);*/

// To Hans: Check the rank of selected D is the same as the rank of H
  RankComputationOfCirMatrix(D, CheckNodeNum, CheckNodeNum, di,swaprow);
 // RankComputationOfCirMatrixH(H, CheckNodeNum,  BitNodeNum, di,swaprow);
// Directly use D to do gaussian elimination is identical to use last part of H

 /* out1=fopen("di.txt","w");
  for(i=0 ; i<CheckNodeNum ; i++)
  {
	if(i%CirculantSize == 0)
	  fprintf(out1, "\n"); 
    
	fprintf(out1, "\n"); 
    for(j=0 ; j<BitNodeNum ; j++)
	{
	  if(j%CirculantSize == 0)
	    fprintf(out1, " "); 
      //Hans modify: delete the failed column and row to simplify the backward substitution
     // if(di[j]!=0 && D[i*CheckNodeNum + j]==1)D[i*CheckNodeNum + j]=0;
	//  if(di[i]!=0 && D[i*CheckNodeNum + j]==1)D[i*CheckNodeNum + j]=0;
      fprintf(out1, "%d", H[i*BitNodeNum + j]); 
	}
  }
  
  for(i=0 ; i<Weight*CirculantSize ; i++)
  {
    fprintf(out1, "%d ", di[i]); 
  }
  fprintf(out1,"\n");
  for(i=0 ; i<CirculantSize ; i++)
  {
    fprintf(out1, "%d ", swaprow[i]); 
  }
  fclose(out1);*/
 /*
  for(i=0 ; i<CheckNodeNum ; i++)
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  if(QC440S->CN[i].BitNodeIdx[j] >= BitNodeNum-CheckNodeNum)
		D[i*CheckNodeNum + QC440S->CN[i].BitNodeIdx[j]-(BitNodeNum-CheckNodeNum)] = 1 ;
	}
 
*/
// To Hans: Solving equation (5). 
//          AngD is the D matrix plus a column of MiuT and do the gaussian elimination
//          Since the row interchange in gaussian elimination, we interchange row back 
//	        to the origin after BackwardSubstitution
//          In order to simplify BackwardSubstitution, we postpone the row interchange 
//          after BackwardSubstitution
  for(w=0 ; w<CirculantSize ; w++)
  {
  for(i=0 ; i<t-c ; i++)
  {
	for(k=0 ; k<CheckNodeNum ; k++)
	{

	  ZiT[k]=0 ;
	}
	for(r=0 ; r<CheckNodeNum ; r++)
    for(v=0 ; v<CheckNodeNum+1 ; v++)
	{
	  AugD[r*(CheckNodeNum+1)+v] = 0 ;
    }
	for(r=0 ; r<CheckNodeNum ; r++)
      for(v=0 ; v<CheckNodeDeg ; v++)
	  {
	  if(QC440S->CN[r].BitNodeIdx[v] >= BitNodeNum-CheckNodeNum)
		AugD[r*(CheckNodeNum+1) + QC440S->CN[r].BitNodeIdx[v]-(BitNodeNum-CheckNodeNum)] = 1 ;
	  } 
	for(j=0 ; j<c ; j++)
	{
	  AugD[ (j * CirculantSize + (CirculantSize-QC440S->HqcLeader[j][i]+w) % CirculantSize +1 )*(CheckNodeNum+1)-1 ] = 1 ;
      //Pack MiuT into last column of AugD 
	  //the above zero can obtain the ziT in the paper, but if we modify the zero to 1 or 2 to obtain the second row of Gi.
	  //We find that it is not circulant relationship
	}
	 
	for(k=0 ; k<CirculantSize ; k++)swaprow[k]=0;
	for(k=0 ; k<CheckNodeNum ; k++)diH[k]=0;
      RankComputationOfCirMatrix(AugD, CheckNodeNum, CheckNodeNum+1, diH,swaprow);
      BackwardSubstitution(di,AugD,ZiT,D,swaprow);


	for(r=0 ; r<CheckNodeNum ; r++)
	{
      if(di[r]==0)
	  {
         GLeadingvector[w*CheckNodeNum*(t-c)+r+i*CheckNodeNum]=ZiT[r];
	  }
	}

  } // End of for(i=0 ; i<t-c ; i++)
  }// End of for(l=0 ; l<CirculantSize ; l++)

 /* out1=fopen("Gqc.txt","a+t");
  fprintf(out1,"\n");
  for(l=0 ; l<CirculantSize ; l++)
  {
  for(i=0 ; i<t-c ; i++)
  {  
	  r=0;
	fprintf(out1, "\n"); 
    for(j=0 ; j<CheckNodeNum ; j++)
	{
      fprintf(out1, "%d ", GLeadingvector[l*CheckNodeNum*(t-c)+j+i*CheckNodeNum]); 
	  if(GLeadingvector[l*CheckNodeNum*(t-c)+j+i*CheckNodeNum]==1)r++;
	}
	fprintf(out1,"(%d)",r);
  }
  	fprintf(out1, "\n"); 
  }
  	fprintf(out1, "\n"); 
  fclose(out1);*/
 

  
  //Hans modify: test the HGT?=0

 /* for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
  {
	  testGH=0;
    for(j=0 ; j<t-c ; j++)
	{	 
		testGH=0;
      for(k=t-c ; k<BitNodeNum/CirculantSize ; k++)
	  {
		  if(GLeadingvector[(QC440S->HqcLeader[i][k]+0)%CirculantSize+(k-t+c)*CirculantSize+j*CheckNodeNum]==1) testGH^=1 ;
		 
	  }
       if(QC440S->HqcLeader[i][j]==0) testGH^=1 ;
	   if(testGH==1)
	   {
		   printf("fail%d\n",i);
	   }
    }

  }// End of for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
  */
 // system("pause");
}




int ParityCheck(Ldpc *QC436)
{
  int i, j;
  int Temp;
  int BitNodeIdx;
  int CheckResult=0;

  for(i=0 ; i<CheckNodeNum && CheckResult==0 ; i++)
  {
	Temp=0 ;
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx =QC436->CN[i].BitNodeIdx[j] ;
	  Temp ^= QC436->BN[BitNodeIdx].TempHD ;
	}
    CheckResult = Temp ;  // If CheckResult = Temp = 1 => Invalid Codeword
  }
  return CheckResult ; // 1 => Invalid,  0 => Valid
}


void CheckNodeProcessing(int layerindex, Ldpc *QC436)
{
  int i,j;
  int BitNodeIdx;
  int ExtArrayIdx;
  double Mji, AbsMji;
  double MinMji, SecMinMji;
  int MinIndex;
  int TotalSign;

  for(i=layerindex*CirculantSize ; i<(layerindex+1)*CirculantSize ; i++)/***************/
  {
	TotalSign = 1 ;
	
	//Min-Sum
	MinMji = Inf;
	MinIndex = 0 ;
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx = QC436->CN[i].BitNodeIdx[j] ;
	  ExtArrayIdx = QC436->BN[BitNodeIdx].index ;
	  Mji = QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx] ;
	  AbsMji = (Mji > 0) ? Mji : -Mji ;
	  if(AbsMji < MinMji)
	  {
		SecMinMji = MinMji;
		MinMji = AbsMji;
		MinIndex = j ;
	  }
	  else if(AbsMji < SecMinMji)
		SecMinMji = AbsMji;

	  TotalSign *= SignFunction(Mji) ;
    }

    // exclude jth branch and compute the extrinsic information
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx = QC436->CN[i].BitNodeIdx[j] ;
	  ExtArrayIdx = QC436->BN[BitNodeIdx].index ;

	  if(j == MinIndex)
	    QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx] = (TotalSign * SignFunction(QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx])) * OffsetFactor  * SecMinMji ;
	  else
		QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx] = (TotalSign * SignFunction(QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx])) * OffsetFactor  * MinMji ;

	  QC436->BN[BitNodeIdx].index = (QC436->BN[BitNodeIdx].index+1) % BitNodeDeg ;
	}
  }
}

void BitNodeProcessing(int layerindex,Ldpc *QC436)
{
  int i,j ;
  double TotalLlr ;
// Update the bit node related to the next check node process .
  for(i=0 ; i<BitNodeNum ; i++)
  {
	TotalLlr = QC436->BN[i].RecLlr ;
	for(j=0 ; j<BitNodeDeg; j++) 
	  TotalLlr += QC436->BN[i].ExtInf[j];
    
	 if( layerindex==(CheckNodeNum/LayerSize-1) ) QC436->BN[i].TempHD = (TotalLlr < 0) ? 1 : 0 ;

	// Update the bit node related to the next check node process .
          j=(layerindex+1) % BitNodeDeg ;  
	  QC436->BN[i].ExtInf[j] = TotalLlr - QC436->BN[i].ExtInf[j];
          
  }
}

void main(void)
{
	  
  int i,j,k,temp,CheckResult,ci;	
  int t=BitNodeNum/CirculantSize, c=CheckNodeNum/CirculantSize;
  Ldpc QC436;
  long seed = -1;
  FILE *out1;
  double SNR;         /* in dB*/
  double var;
  int IterNum;
  int BlockNum;
  int ErrorFlag, ErrorBlockNumber, ErrorBits, RawErrorBits;
  int ParityCheckResult;


  RandomGenHMatrix(&QC436) ;
  ComputeEncoder(&QC436,GLeadingvector) ; 

  out1=fopen("FinalResult.txt","a+t");
  fprintf(out1, "\nQC-LDPC with %d Bit Nodes (including %d Shorten Bit) and %d Check Nodes and Circulant Size = %d", BitNodeNum, ShortenNum, CheckNodeNum, CirculantSize);
  fprintf(out1, "\nOffset Min-Sum Decoding Algorithm with factor %f and double input and Maximum %d Iterations", OffsetFactor, MaxIterNum);
  fprintf(out1, "\nBit Node Degree = %d, Check Node Degree = %d, with Random H Matrix with Seed = %d", BitNodeDeg, CheckNodeDeg, InitialMatrixSeed);
  fprintf(out1, "\n\n") ;
  fclose(out1);

  //Random generate data
  	seed = 131073;
	srand(seed);
  	for (i = 0; i < BitNodeNum-CheckNodeNum; i++){data[i] =0;
	//	data[i] =  rand() % 2;
	}
  //Encode data to parity bb 

		
    for (j = 0; j < CheckNodeNum; j++)
	{
		temp=0;
	  for (k = 0; k < BitNodeNum-CheckNodeNum; k++)
	  {
		  if(GLeadingvector[(k % CirculantSize)*CheckNodeNum*(t-c)+(int)(k / CirculantSize)*CheckNodeNum+j]==1 && data[k]==1) temp^=1;

	  }
		bb[j]=temp;
        recd[j + BitNodeNum-CheckNodeNum] =  bb[j]; 
	}
	for (i = 0; i < BitNodeNum-CheckNodeNum; i++)
		recd[i] = data[i];

/* //Check parity correct 
  for(i=0 ; i<CheckNodeNum && CheckResult==0 ; i++)
  {
	temp=0 ;
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx = QC436.CN[i].BitNodeIdx[j] ;
	  temp ^= recd[BitNodeIdx] ;
	}
    CheckResult = temp ; 
	// If CheckResult = Temp = 1 => Invalid Codeword
    if(temp==1)break;
  }
*/
  //Assume All Zero Codeword Is Transmitted so that we can ignore LDPC encoder at this stage

	for(SNR=3 ; SNR<=4.1 ; SNR+=0.1)
  {
	printf("\nSNR=%f", SNR) ;
  //  var=(double)(1.0/2.0*pow(10,-(double)SNR/10));
      var=(double)((BitNodeNum)/(BitNodeNum-CheckNodeNum)/2.0*pow(10,-(double)SNR/10));
    ErrorBlockNumber = 0 ;
	ErrorBits = 0 ;
	RawErrorBits = 0 ;
	//for(BlockNum=0 ; ErrorBlockNumber<50 ; BlockNum++)
	for(BlockNum=0 ; ErrorBlockNumber<50; BlockNum++)
	{
	  if(BlockNum % 10000 == 0)
	  {
	    printf(".") ;
		out1=fopen("Result.txt","w");
        fprintf(out1, "\nSNR:%f, RawBER=%d/%d=%e, BER=%d/%d=%e, BLER=%d/%d=%e\n" , SNR, RawErrorBits, (BitNodeNum-ShortenNum)*BlockNum, (double)RawErrorBits/((BitNodeNum-ShortenNum)*BlockNum), 
		  ErrorBits, (BitNodeNum-ShortenNum)*BlockNum, (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum), ErrorBlockNumber, BlockNum ,(double)ErrorBlockNumber/(BlockNum) );
        fclose(out1);
	  }


	  // Channel Noise and Initialization
	  
	  // Part 1: Shorten bits are not transmitted and are known to the receiver
	  for(i=0 ; i<ShortenNum ; i++)
	  {
		//
        QC436.BN[i].RecLlr = (double)(2 * (Inf) / var) ;
		QC436.BN[i].TempHD = (QC436.BN[i].RecLlr <0 ) ? 1 : 0 ;
	    for(j=0 ; j<BitNodeDeg ; j++)
	      QC436.BN[i].ExtInf[j] = QC436.BN[i].RecLlr ;
		QC436.BN[i].index = 0 ;     
	  }

	  // Part 2: Real transmitted bits
	  for(i=ShortenNum ; i<BitNodeNum ; i++)
	  {
      //  QC436.BN[i].RecLlr = (double)(2 * ( (1-2*recd[i]) + sqrt( var)*gaussian_noise(&seed) ) / var) ;
		  QC436.BN[i].RecLlr = (double)(2 * ((1-2*recd[i])) / var) ;
		QC436.BN[i].TempHD = (QC436.BN[i].RecLlr <0 ) ? 1 : 0 ;
	    for(j=0 ; j<BitNodeDeg ; j++)
	      QC436.BN[i].ExtInf[j] = QC436.BN[i].RecLlr ;
		QC436.BN[i].index = 0 ;     
		
		if(QC436.BN[i].TempHD != recd[i])
	      RawErrorBits ++ ;
                  for(j=0 ; j<BitNodeDeg ; j++) QC436.BN[i].ExtInf[j] =0 ;// For 1st iteration layer decoding	
	  }

	  IterNum = 0 ;
	  ParityCheckResult = ParityCheck(&QC436) ;
	  for(IterNum=0 ; IterNum<MaxIterNum && ParityCheckResult==1 ; IterNum++)
         {    
        for( ci=0;ci<(CheckNodeNum/LayerSize);ci++)/***************/
		{
               /*Since 1st iteration next layer has no update ExtInf, 
                 it can not load all channel value into ExtInf at the initial step. 
                 BER will result in only 0.1 difference. */ 
                 if(IterNum==0)for(i=0 ; i<BitNodeNum ; i++)QC436.BN[i].ExtInf[ci] = QC436.BN[i].RecLlr ;
                  CheckNodeProcessing(ci,&QC436) ;

           BitNodeProcessing(ci,&QC436) ;
		}
	   ParityCheckResult = ParityCheck(&QC436) ;
      } // End of for(IterNum=0 ; IterNum<MaxIterNum ; IterNum++)
	  ErrorFlag=0;
      for(i=ShortenNum ; i<BitNodeNum ; i++)
	  {
        if(QC436.BN[i].TempHD != recd[i])
		{
	      ErrorBits++;
          ErrorFlag=1;
		}
	  }
      if(ErrorFlag==1)
	    ErrorBlockNumber++;

      /*Report the error pattern as trapping set*/
	  if(ErrorFlag==1)
	  {
        
        out1=fopen("TrappingsetResult.txt","a+t");
        fprintf(out1,"(%d)th ",ErrorBlockNumber); 
		for(i=ShortenNum ; i<BitNodeNum ; i++)
		{
          if(QC436.BN[i].TempHD^recd[i]!=0)
	         fprintf(out1,"%d  ",i); 
	   
		}
        fprintf(out1,"\n");
		fclose(out1);
	  } 
    } // End of for(BlockNum=0 ; BlockNum<100 ; BlockNum++)

    out1=fopen("FinalResult.txt","a+t");
    fprintf(out1, "\nSNR:%f, var=%f, RawBER=%d/%d=%e, BER=%d/%d=%e, BLER=%d/%d=%e\n" , SNR, var, RawErrorBits, (BitNodeNum-ShortenNum)*BlockNum, (double)RawErrorBits/((BitNodeNum-ShortenNum)*BlockNum), 
		ErrorBits, (BitNodeNum-ShortenNum)*BlockNum, (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum), ErrorBlockNumber, BlockNum ,(double)ErrorBlockNumber/(BlockNum) );
    fclose(out1);

    out1=fopen("BER.txt","a+t");
    fprintf(out1, "%e ", (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum));
    fclose(out1);
    printf( "%e ", (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum));
    out1=fopen("BLER.txt","a+t");
    fprintf(out1, "%e ", (double)ErrorBlockNumber/(BlockNum));
    fclose(out1);

  } // End of for(SNR=1.0 ; SNR<=5.0 ; SNR+=1.0)

}
