#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <time.h>

#define CirculantSize 256
#define Weight 4
#define InitialMatrixSeed 123

#define RateParameter 10 // 8 code rate = 1 - (1/RateParameter) ;

#define BitNodeNum (CirculantSize*RateParameter*Weight)  
#define CheckNodeNum (CirculantSize*Weight)
#define ShortenNum 0

#define OffsetFactor 0.75

#define BitNodeDeg Weight
#define CheckNodeDeg (Weight*RateParameter)
#define MaxIterNum 16

#define Inf 1000000
int recd[1048576],recd1[1048576], data[1048576], bb[548576];
int seed;
int temparray[CheckNodeNum];
int Gen[BitNodeNum-CheckNodeNum][CheckNodeNum];
int GLeadingvector[CheckNodeDeg-BitNodeDeg][CheckNodeNum];
int H[CheckNodeNum*BitNodeNum];
int Hpi[CheckNodeNum][BitNodeNum];
int D[CheckNodeNum*CheckNodeNum];
int AugD[CheckNodeNum*(CheckNodeNum+1)];
int AugH[CheckNodeNum*(BitNodeNum+1)];
int di[CheckNodeNum+1]={0};
int diH[CheckNodeNum+1]={0};
int swaprow[50*CheckNodeNum]={0};int swaprowH[50*CheckNodeNum]={0};
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

void RandomGenHMatrix(Ldpc *QC440S)
{
  int i,j;
  int k;
  long MatrixSeed = InitialMatrixSeed;
  FILE *out1;

  for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{
	  if(i!=0 && j!=0) // Not in the first column or the first row
	    QC440S->HqcLeader[i][j] = (int) (ran0(&MatrixSeed) * CirculantSize) ;
 	  else
		QC440S->HqcLeader[i][j] = 0 ;

	  for(k=0 ; k<CirculantSize ; k++)
	    QC440S->CN[i*CirculantSize+k].BitNodeIdx[j] = j*CirculantSize+ ( (QC440S->HqcLeader[i][j]+k)%CirculantSize ) ;
    }
/*	j=0;
	
    for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
	for(k = 0; k < BitNodeNum/CirculantSize; k++)
	{

		QC440S->HqcLeader2[i][k] = QC440S->HqcLeader[i][j];
	   j++;
	   	if(j==0 || j==5 || j==11 ||j==15)j++;
     
	}
    for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
	{
     	QC440S->HqcLeader2[i][0] = QC440S->HqcLeader[i][(CheckNodeDeg-4)];

     	QC440S->HqcLeader2[i][5] = QC440S->HqcLeader[i][(CheckNodeDeg-3)];

     	QC440S->HqcLeader2[i][11] = QC440S->HqcLeader[i][(CheckNodeDeg-2)];

     	QC440S->HqcLeader2[i][15] = QC440S->HqcLeader[i][(CheckNodeDeg-1)];
	}*/

  out1=fopen("RandomHMatrix.txt","a+t");
  for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
  {
    fprintf(out1, "\n"); 
    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{
       fprintf(out1, "%3d ", QC440S->HqcLeader[i][j]); 
	}
  }
  fclose(out1);
  

}

int ParityCheck(Ldpc *QC440S)
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
	  BitNodeIdx = QC440S->CN[i].BitNodeIdx[j] ;
	  Temp ^= QC440S->BN[BitNodeIdx].TempHD ;
	}
    CheckResult = Temp ;  // If CheckResult = Temp = 1 => Invalid Codeword
  }
  return CheckResult ; // 1 => Invalid,  0 => Valid
}

void CheckNodeProcessing(Ldpc *QC440S)
{
  int i,j;
  int BitNodeIdx;
  int ExtArrayIdx;
  double Mji, AbsMji;
  double MinMji, SecMinMji;
  int MinIndex;
  int TotalSign;

  for(i=0 ; i<CheckNodeNum ; i++)
  {
	TotalSign = 1 ;
	
	//Min-Sum
	MinMji = Inf;
	MinIndex = 0 ;
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx = QC440S->CN[i].BitNodeIdx[j] ;
	  ExtArrayIdx = QC440S->BN[BitNodeIdx].index ;
	  Mji = QC440S->BN[BitNodeIdx].ExtInf[ExtArrayIdx] ;
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
	  BitNodeIdx = QC440S->CN[i].BitNodeIdx[j] ;
	  ExtArrayIdx = QC440S->BN[BitNodeIdx].index ;

	  if(j == MinIndex)
	    QC440S->BN[BitNodeIdx].ExtInf[ExtArrayIdx] = (TotalSign * SignFunction(QC440S->BN[BitNodeIdx].ExtInf[ExtArrayIdx])) * OffsetFactor * SecMinMji ;
	  else
		QC440S->BN[BitNodeIdx].ExtInf[ExtArrayIdx] = (TotalSign * SignFunction(QC440S->BN[BitNodeIdx].ExtInf[ExtArrayIdx])) * OffsetFactor * MinMji ;

	  QC440S->BN[BitNodeIdx].index = (QC440S->BN[BitNodeIdx].index+1) % BitNodeDeg ;
	}
  }
}

void BitNodeProcessing(Ldpc *QC440S)
{
  int i,j ;
  double TotalLlr ;
  for(i=0 ; i<BitNodeNum ; i++)
  {
	TotalLlr = QC440S->BN[i].RecLlr ;
	for(j=0 ; j<BitNodeDeg ; j++)
	  TotalLlr += QC440S->BN[i].ExtInf[j];
    
	QC440S->BN[i].TempHD = (TotalLlr < 0) ? 1 : 0 ;

	for(j=0 ; j<BitNodeDeg ; j++)
	  QC440S->BN[i].ExtInf[j] = TotalLlr - QC440S->BN[i].ExtInf[j];
  }
}

void main(void)
{
	  
  int i,j,k,u,temp,CheckResult,BitNodeIdx;	
  int t=BitNodeNum/CirculantSize, c=CheckNodeNum/CirculantSize;
  Ldpc QC440S;
  long seed = -1;
  
 
  char buffer[10000];
  double SNR;         /* in dB*/
  double var;
  int IterNum;
  int BlockNum;
  int ErrorFlag, ErrorBlockNumber, ErrorBits, RawErrorBits;
  int ParityCheckResult;
  FILE *out1;
  FILE *out2;
	out2=fopen("Generater.txt","r");    
	fgets(buffer,10000,out2);
	
	//	Load data_cv
	for (i=0 ; i<CheckNodeDeg-BitNodeDeg ; i++){
		for (j=0 ; j<CheckNodeNum ; j++){
			GLeadingvector[i][j]=buffer[j]-48;
		}
		fgets(buffer,10000,out2);	
	}
	fclose(out2);

  RandomGenHMatrix(&QC440S) ;
 // ComputeEncoder(&QC440S,GLeadingvector) ; 
   for(i =0; i <CheckNodeDeg-BitNodeDeg ; i++)
   {
   for (j = 0; j < CheckNodeNum; j++)
   {
	   temparray[j]=GLeadingvector[i][j];
   }
   for (j = 0; j < CheckNodeNum-1; j++)
   {
		Gen[i*CirculantSize][j]=temparray[j];
   }
   for (k = 1; k < CirculantSize; k++)	  
   {
       Gen[k+i*CirculantSize][0]=temparray[CirculantSize-1];
	   Gen[k+i*CirculantSize][CirculantSize]=temparray[2*CirculantSize-1];
	   Gen[k+i*CirculantSize][2*CirculantSize]=temparray[3*CirculantSize-1];
	   Gen[k+i*CirculantSize][3*CirculantSize]=temparray[4*CirculantSize-1];
      for ( u =0; u< Weight; u++)
	  for (j = 0; j < CirculantSize-1; j++)
	  {
		  Gen[k+i*CirculantSize][j+1+u*CirculantSize]=temparray[j+u*CirculantSize];
	  }
	   for (j = 0; j < CheckNodeNum; j++)
           temparray[j]=Gen[k+i*CirculantSize][j];
   }
   }
	


  //Random generate data
  	seed = 131073;
	srand(seed);
  for(SNR=3 ; SNR<=4 ; SNR+=0.1)
  {
  	for (i = 0; i < BitNodeNum-CheckNodeNum; i++){data[i] =0;
	//	data[i] =  rand() % 2;
	}
  //Encode data to parity bb 

		
    for (j = 0; j < CheckNodeNum; j++)
	{
		temp=0;
	  for (k = 0; k < BitNodeNum-CheckNodeNum; k++)
	  {
         if(Gen[k][j]==1 && data[k]==1) temp^=1;

	  }
		bb[j]=temp;
        recd[j + BitNodeNum-CheckNodeNum] =  bb[j]; 
	}
	for (i = 0; i < BitNodeNum-CheckNodeNum; i++)
		recd[i] = data[i];

	j=0;
	for(k = 0; k < CheckNodeDeg; k++)
	{
      if(k==1 || k==26 || k==33 ||k==37)k++;
	//	if(k==0 || k==5 || k==11 ||k==15)k++;
	  for (i = 0; i < CirculantSize; i++)
		recd1[k*CirculantSize+i] = recd[j*CirculantSize+i];
	   j++;
     
	}

    for (i = 0; i < CirculantSize; i++)
     	recd1[1*CirculantSize+i] = recd[(CheckNodeDeg-4)*CirculantSize+i];
    for (i = 0; i < CirculantSize; i++)
     	recd1[26*CirculantSize+i] = recd[(CheckNodeDeg-3)*CirculantSize+i];
    for (i = 0; i < CirculantSize; i++)
     	recd1[33*CirculantSize+i] = recd[(CheckNodeDeg-2)*CirculantSize+i];
    for (i = 0; i < CirculantSize; i++)
     	recd1[37*CirculantSize+i] = recd[(CheckNodeDeg-1)*CirculantSize+i];

  CheckResult=0;


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
        QC440S.BN[i].RecLlr = (double)(2 * (Inf) / var) ;
		QC440S.BN[i].TempHD = (QC440S.BN[i].RecLlr <0 ) ? 1 : 0 ;
	    for(j=0 ; j<BitNodeDeg ; j++)
	      QC440S.BN[i].ExtInf[j] = QC440S.BN[i].RecLlr ;
		QC440S.BN[i].index = 0 ;     
	  }

	  // Part 2: Real transmitted bits
	  for(i=ShortenNum ; i<BitNodeNum ; i++)
	  {
      //  QC440S.BN[i].RecLlr = (double)(2 * ((1-2*recd1[i]) + sqrt( var)*gaussian_noise(&seed)) ) ;
		  QC440S.BN[i].RecLlr = (double)(2 * ((1-2*recd1[i])) / var) ;
		QC440S.BN[i].TempHD = (QC440S.BN[i].RecLlr <0 ) ? 1 : 0 ;
	    for(j=0 ; j<BitNodeDeg ; j++)
	      QC440S.BN[i].ExtInf[j] = QC440S.BN[i].RecLlr ;
		QC440S.BN[i].index = 0 ;     
		
		if(QC440S.BN[i].TempHD != recd[i])
	      RawErrorBits ++ ;
	  }

	  IterNum = 0 ;
	  ParityCheckResult = ParityCheck(&QC440S) ;
	  for(IterNum=0 ; IterNum<MaxIterNum && ParityCheckResult==1 ; IterNum++)
	  {
        CheckNodeProcessing(&QC440S) ;
        BitNodeProcessing(&QC440S) ;
		ParityCheckResult = ParityCheck(&QC440S) ;
	  } // End of for(IterNum=0 ; IterNum<MaxIterNum ; IterNum++)

      
	  ErrorFlag=0;
      for(i=ShortenNum ; i<BitNodeNum ; i++)
	  {
        if(QC440S.BN[i].TempHD != recd1[i])
		{
	      ErrorBits++;
          ErrorFlag=1;
		  

		}
	  }

      if(ErrorFlag==1)
	    ErrorBlockNumber++;
	 
	  /*Report the error pattern as trapping set*/
	  /*if(ErrorFlag==1)
	  {
        
        out1=fopen("TrappingsetResult.txt","a+t");
        fprintf(out1,"(%d)th  ",ErrorBlockNumber); 
		for(i=ShortenNum ; i<BitNodeNum ; i++)
		{
          if(QC440S.BN[i].TempHD^recd[i]!=0)
	         fprintf(out1,"%d  ",i); 
	   
		}
        fprintf(out1,"\n");
		fclose(out1);
	  } */
 
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