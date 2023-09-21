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
int             recd[1048576],recd1[1048576], data[1048576], bb[548576];
int             seed;
int temparray[CheckNodeNum];
int GLeadingvector[CheckNodeDeg-BitNodeDeg][CheckNodeNum];
int Gen[BitNodeNum-CheckNodeNum][CheckNodeNum];
char H[CheckNodeNum*BitNodeNum];
int Hpi[CheckNodeNum][BitNodeNum];
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
  int HqcLeader2[CheckNodeNum/CirculantSize][BitNodeNum/CirculantSize];
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

//	  for(k=0 ; k<CirculantSize ; k++)
//	    QC440S->CN[i*CirculantSize+k].BitNodeIdx[j] = j*CirculantSize+ ( (QC440S->HqcLeader[i][j]+k)%CirculantSize ) ;
    }
	j=0;
	

	for(k = 0; k < CheckNodeDeg-4; k++)
	{
       for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
		QC440S->HqcLeader2[i][k] = QC440S->HqcLeader[i][j];
	   j++;
	   	if(j==1 || j==26 || j==33 ||j==37)j++;
	   
	}
    for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
	{
     	QC440S->HqcLeader2[i][(CheckNodeDeg-4)] = QC440S->HqcLeader[i][1];

     	QC440S->HqcLeader2[i][(CheckNodeDeg-3)] = QC440S->HqcLeader[i][26];

     	QC440S->HqcLeader2[i][(CheckNodeDeg-2)] = QC440S->HqcLeader[i][33];

     	QC440S->HqcLeader2[i][(CheckNodeDeg-1)] = QC440S->HqcLeader[i][37];
	}
  for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
   for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{

	  for(k=0 ; k<CirculantSize ; k++)
	    QC440S->CN[i*CirculantSize+k].BitNodeIdx[j] = j*CirculantSize+ ( (QC440S->HqcLeader2[i][j]+k)%CirculantSize ) ;
    }
  out1=fopen("RandomHMatrix.txt","a+t");
  for(i=0 ; i<CheckNodeNum/CirculantSize ; i++)
  {
    fprintf(out1, "\n"); 
    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{
       fprintf(out1, "%3d ", QC440S->HqcLeader2[i][j]); 
	}
  }
  fclose(out1);
  

}
int swapnum;
	int i,j,r,v,k,l,temp,w;int t, c;char tmp;




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
	  else if(AbsMji < SecMinMji && AbsMji!=MinMji)
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
     /**quantization**/
//	TotalLlr=(TotalLlr>3)?3:TotalLlr;
//	TotalLlr=(TotalLlr<-3)?-3:TotalLlr;
	 /**quantization**/
	 if( layerindex==(CheckNodeNum/LayerSize-1) ) QC436->BN[i].TempHD = (TotalLlr < 0) ? 1 : 0 ;

	// Update the bit node related to the next check node process .
          j=(layerindex+1) % BitNodeDeg ;  
	  QC436->BN[i].ExtInf[j] = TotalLlr - QC436->BN[i].ExtInf[j];
	  /**quantization**/
//	  QC436->BN[i].ExtInf[j]=(QC436->BN[i].ExtInf[j]>10)?10:QC436->BN[i].ExtInf[j];
//	  QC436->BN[i].ExtInf[j]=(QC436->BN[i].ExtInf[j]<-10)?-10:QC436->BN[i].ExtInf[j];
	  /**quantization**/
          
  }
}

void main(void)
{
	  
  int i,j,k,u,temp,CheckResult,ci;	
  int t=BitNodeNum/CirculantSize, c=CheckNodeNum/CirculantSize;
  Ldpc QC436;
  long seed = -1;
    int FAb[10],VAb[10],FVb[20];
	  char buffer[10000];
  FILE *out1;
  double SNR;         /* in dB*/
  double var,FA,FA_temp,ss,VA,VA_temp,VA_temp1,VA_temp2,FV,FV_temp;
  int Sgn,FA_int,VA_int,FV_int;
  int IterNum;
  int BlockNum;
  int ErrorFlag, ErrorBlockNumber, ErrorBits, RawErrorBits;
  int ParityCheckResult;

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

  RandomGenHMatrix(&QC436) ;
 // ComputeEncoder(&QC440S,GLeadingvector) ; 
  //沒有再做行交換只是把G產生出來 
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
  for(SNR=3 ; SNR<=7; SNR+=0.2)
  {	
	printf("\nSNR=%f", SNR) ;
  //  var=(double)(1.0/2.0*pow(10,-(double)SNR/10));
      var=(double)((BitNodeNum)/(BitNodeNum-CheckNodeNum)/2.0*pow(10,-(double)SNR/10));
          /**quantization**/
	      VA=sqrt(var)/var;
		  VA_temp=VA;
		  VA_int=VA;
          ss=1;
  	    for(j=0 ; j<6 ; j++)
		{
			ss=ss/2;
		  VAb[j]=(VA-ss>0)?1:0;
		  if(VA-ss>0)VA=VA-ss;
		}		 
        //  printf("%f %d%d%d\n",VA_temp,VAb[0],VAb[1],VAb[2]);
        VA_temp1=VA_temp-VA;
        /**quantization**/

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

  	for (i = 0; i < BitNodeNum-CheckNodeNum; i++){//data[i] =0;
		data[i] =  rand() % 2;
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


  CheckResult=0;

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
	out1=fopen("codeword1.txt","a+t");  
	  // Part 2: Real transmitted bits
	  for(i=ShortenNum ; i<BitNodeNum ; i++)
	  {
		  /**quantization**/
		   FA=gaussian_noise(&seed);
		  FA_temp=FA;
		  FA_int=FA;
		  Sgn=(FA<0)?1:0;
		  FA=(FA<0)?FA*(-1):FA;
          ss=1;
  	    for(j=0 ; j<6 ; j++)
		{
			ss=ss/2;
		  if(Sgn==0) FAb[j]=(FA-ss>=0)?1:0;
          else FAb[j]=(FA-ss>=0)?0:1; 
		  if(FA-ss>=0)FA=FA-ss;
		}	
		FA_temp=(FA_temp<0)?FA_temp+FA:FA_temp-FA;
		/**quantization**/
		 //  printf("%f %d%d%d\n",FA_temp,FAb[0],FAb[1],FAb[2]);
		//  QC436.BN[i].RecLlr = (double)(2 * ((1-2*recd1[i]) + sqrt( var)*gaussian_noise(&seed)) ) ;
	//	  QC436.BN[i].RecLlr = (double)(  (1-2*recd[i])+ VA*gaussian_noise(&seed)  ) ;
	     QC436.BN[i].RecLlr = (double)(2 * ( (1-2*recd[i]) + sqrt( var)*gaussian_noise(&seed) ) / var) ;  
		/**quantization**/  
		  FA=QC436.BN[i].RecLlr;
		  FA_temp=FA;
		  FA_int=FA;
		  Sgn=(FA<0)?1:0;
		  FA=(FA<0)?FA*(-1):FA;
          ss=1;
  	    for(j=0 ; j<2 ; j++)
		{
			ss=ss/2;
		  if(Sgn==0) FAb[j]=(FA-ss>=0)?1:0;
          else FAb[j]=(FA-ss>=0)?0:1; 
		  if(FA-ss>=0)FA=FA-ss;
		}	
		FA_temp=(FA_temp<0)?FA_temp+FA:FA_temp-FA;
		QC436.BN[i].RecLlr =FA_temp;
		/**quantization**/
	//	if(i%16==0)fprintf(out1,"\n");/*********************************************************************************/
   //    fprintf(out1,"%d%d%d",Sgn,FAb[0],FAb[1]);/*********************************************************************************/
		 // QC436.BN[i].RecLlr = (double)(2 * ((1-2*recd[i])) / var) ;
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
	//	if(i%16==0)fprintf(out1,"\n");/*********************************************************************************/
	  //  fprintf(out1,"%d11",recd[i]);/*********************************************************************************/

        if(QC436.BN[i].TempHD != recd[i])
		{
		/**quantization**/  
		  FA=i;
          ss=8192;
  	    for(j=0 ; j<14 ; j++)
		{			
			FAb[j]=(FA-ss>=0)?1:0; 
		  if(FA-ss>=0)FA=FA-ss;
		  ss=ss/2;
		//  fprintf(out1,"%d",FAb[j]);/*********************************************************************************/
		}	
    //	fprintf(out1,"\n");/*********************************************************************************/
		/**quantization**/
	      ErrorBits++;
          ErrorFlag=1;
		}
	  }
      if(ErrorFlag==1)
	    ErrorBlockNumber++;
fclose(out1);
      /*Report the error pattern as trapping set*/

    } // End of for(BlockNum=0 ; BlockNum<100 ; BlockNum++)

    out1=fopen("FinalResult.txt","a+t");
    fprintf(out1, "\nSNR:%f, var=%f, RawBER=%d/%d=%e, BER=%d/%d=%e, BLER=%d/%d=%e\n" , SNR, var, RawErrorBits, (BitNodeNum-ShortenNum)*BlockNum, (double)RawErrorBits/((BitNodeNum-ShortenNum)*BlockNum), 
		ErrorBits, (BitNodeNum-ShortenNum)*BlockNum, (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum), ErrorBlockNumber, BlockNum ,(double)ErrorBlockNumber/(BlockNum) );
    fclose(out1);

    out1=fopen("BER.txt","a+t");
    fprintf(out1, "%e ", (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum));
    fclose(out1);
    printf( "%e   %e ",(double)RawErrorBits/((BitNodeNum-ShortenNum)*BlockNum), (double)ErrorBits/((BitNodeNum-ShortenNum)*BlockNum));
    out1=fopen("BLER.txt","a+t");
    fprintf(out1, "%e ", (double)ErrorBlockNumber/(BlockNum));
    fclose(out1);

  } // End of for(SNR=1.0 ; SNR<=5.0 ; SNR+=1.0)

}
