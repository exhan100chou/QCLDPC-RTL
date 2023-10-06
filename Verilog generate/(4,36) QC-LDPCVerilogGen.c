#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BitNodeNum 10240 //9216
#define CheckNodeNum 1024 //1024
#define CirculantSize 256 //256

#define InitialMatrixSeed 123

#define BitNodeDeg 4
#define CheckNodeDeg 40
#define MaxIterNum 16

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
}Ldpc;
#define tabnum 3
typedef struct
{
   int col_indx[CheckNodeDeg];
   int row_indx[CheckNodeDeg];

}Tabindx;
typedef struct
{
  Tabindx TI[tabnum];

}QCTable;
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

int cmn_comp(int tabindx,QCTable *tab36)
{
    int i,j,k,cmn;
   for(i=0,cmn=0 ; i<BitNodeNum/CirculantSize && cmn==0 ; i++)
   {

	for(j=0; j<BitNodeNum/CirculantSize ; j++)
	{

		for(k=0 ; k<BitNodeNum/CirculantSize ; k++)
		{
        if(k!=i&&j!=i)
          if( (tab36->TI[tabindx].col_indx[i]+tab36->TI[tabindx].row_indx[j])%CirculantSize==(tab36->TI[tabindx].col_indx[k]+tab36->TI[tabindx].row_indx[i])%CirculantSize )
		  {
			cmn++;
		  }
		}
	}
	if(cmn!=0)printf("common value in %dth column and row!!\n",i);
	
   }
 return(cmn);
}
QCTable tab36;
void FFGenHMatrix(Ldpc *QC436)
{
  int i,j,t=0,ir=0,fr;
  int k,cmn;

  int CirculantLeader[CheckNodeNum/CirculantSize][BitNodeNum/CirculantSize]={0} ;
  long MatrixSeed = InitialMatrixSeed;
  int Hnum=1;
  FILE *out1;

  for(i=ir ; i<CheckNodeNum/CirculantSize ; i++)
    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
	{
	  if(i!=0 && j!=0) // Not in the first column or the first row
	  {
			//  CirculantLeader[i][j] = (int)(((int)( pow(3.0,(double)i) * pow(5.0,(double)j)) )%CirculantSize) ;
		  CirculantLeader[i][j] = (int) (ran0(&MatrixSeed) * CirculantSize) ;

	  }
      
	  for(k=0 ; k<CirculantSize ; k++)
	    QC436->CN[i*CirculantSize+k].BitNodeIdx[j] = j*CirculantSize+ ( (CirculantLeader[i][j]+k)%CirculantSize ) ;
    }
 
  out1=fopen("RandomHMatrix.txt","w");
  cmn=0;
for(fr=1 ; fr<=5; fr++)
{
  fprintf(out1, "\n");
  for(i=1 ; i<CheckNodeNum/CirculantSize ; i++)
  {
 
//	s=4;
//	ir=0;
    for(j=8*(fr-1) ; j<8*fr; j++)
	{
     //  fprintf(out1, "%3d ", CirculantLeader[i][j]); 
	 //  fprintf(out1, "%d      %d       %d\n",CirculantLeader[i][j],  CirculantLeader[i][j]%16,  CirculantLeader[i][j]/16);
	 fprintf(out1, "%d       %d\n", CirculantLeader[i][j]%16,  (16-CirculantLeader[i][j]/16)*15);

	 /*   if( CirculantLeader[i][j]/32!=0)
		  fprintf(out1, "rtm%d[%d:0],rtm%d[119:%d]\n",j+1,15*(int)(CirculantLeader[i][j]/32)-1,j+1,15*(int)(CirculantLeader[i][j]/32));
	   else
          fprintf(out1, "rtm%d\n",j);*/
/*
	   if( CirculantLeader[i][j]/32!=0)
		  fprintf(out1, "mtv%d[%d:0],mtv%d[119:%d]\n",j+1,119-15*(int)(CirculantLeader[i][j]/32)-1,j+1,119-15*(int)(CirculantLeader[i][j]/32));
	   else
          fprintf(out1, "mtv%d\n",j);*/
/*
		if(j%10==0&&j!=0&&s>0)
		{
			j=j-10;
		    s--;
			if(s==0)
			{ 
				s=4;
				ir++;
				j=10*ir;
				cmn++;
                fprintf(out1,"(i=%d)(%d)\n",i,cmn);
				if(ir==4)break;
			}
		}
	   sh=CirculantLeader[i][j]/64;

        for(k=4 ; k>0 ; k--)fprintf(out1, "mtv%d[%d:%d]\n",  j%10+1, 5*k-1+12*((4-s+sh)%4),5*(k-1)+12*((4-s+sh)%4));
	*/	
	}
	
  }
}

 /*    for(j=0 ; j<BitNodeNum/CirculantSize ; j++)
   {
	   tab36.TI[0].col_indx[j] =CirculantLeader[1][j];
       tab36.TI[0].row_indx[j] =CirculantLeader[2][j];
	   tab36.TI[1].col_indx[j] =CirculantLeader[1][j];
	   tab36.TI[1].row_indx[j] =CirculantLeader[3][j];
	   tab36.TI[2].col_indx[j] =CirculantLeader[2][j];
	   tab36.TI[2].row_indx[j] =CirculantLeader[3][j];
   }

  for(i=0,cmn=0;i<tabnum&&cmn==0;i++)cmn=cmn_comp(i, &tab36);


   if(cmn!=0)printf("try %dth H\n",t);

 
*/
   for(i=1;i<41;i++)
   {
/*	 if(i%40!=0) fprintf(out1,"\nVNP VNP%d_%d(\n",i%40,i/40+1); else fprintf(out1,"\nVNP VNP%d_%d(\n",40,i/40);
	 fprintf(out1,"         .CLK(CLK),\n");
	 fprintf(out1,"         .RESET_N(RESET_N),\n");
	 fprintf(out1,"         .dec_en(dec_en),\n");
	 fprintf(out1,"         .hb_en(hb_en),\n");
     if(i%40!=0)fprintf(out1,"         .mtv(mtv%d_%d_w),\n",i%40,i/40+1); else fprintf(out1,"         .mtv(mtv%d_%d_w),\n",40,i/40);
	 if(i%40!=0)fprintf(out1,"         .summtv(summtv%d_%d_w),\n",i%40,i/40+1); else fprintf(out1,"         .summtv(summtv%d_%d_w),\n",40,i/40);
	 if(i%40!=0)fprintf(out1,"         .vtc(vtc%d_%d),\n",i%40,i/40+1); else fprintf(out1,"         .vtc(vtc%d_%d),\n",40,i/40); 
	 if(i%40!=0)fprintf(out1,"         .hardbit(hardbit%d_%d)\n",i%40,i/40+1); else fprintf(out1,"         .hardbit(hardbit%d_%d)\n",40,i/40);
	 fprintf(out1,"         );\n ");
            */
     // fprintf(out1,"            mtv%d_%d_w,mtv%d_%d_w,mtv%d_%d_w,mtv%d_%d_w,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1,(i+3)%40,(i+3)/40+1);
     //   fprintf(out1,"            summtv%d_%d_w,summtv%d_%d_w,summtv%d_%d_w,summtv%d_%d_w,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1,(i+3)%40,(i+3)/40+1); 
	  //  fprintf(out1,"            ctm%d_%d,ctm%d_%d,ctm%d_%d,ctm%d_%d,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1,(i+3)%40,(i+3)/40+1);
	//   fprintf(out1,"            vtc%d_%d,vtc%d_%d,vtc%d_%d,vtc%d_%d,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1,(i+3)%40,(i+3)/40+1);
	//   fprintf(out1,"            mtv%d_w,mtv%d_w,mtv%d_w,mtv%d_w,\n",i,i+1,i+2,i+3);
	  //   fprintf(out1,"            rtm%d_%d_w,rtm%d_%d_w,rtm%d_%d_w,rtm%d_%d_w,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1,(i+3)%40,(i+3)/40+1);
     //    fprintf(out1,"            vtm%d_%d,vtm%d_%d,vtm%d_%d,vtm%d_%d,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1,(i+3)%40,(i+3)/40+1); 
	 //  fprintf(out1,"                hardbit%d_%d,hardbit%d_%d,hardbit%d_%d,\n",i%40,i/40+1,(i+1)%40,(i+1)/40+1,(i+2)%40,(i+2)/40+1);  
 /*    fprintf(out1,"wire [59:0] rtm%d=\n",i);                                                                                                                                                                                                                                                                                                       
     fprintf(out1,"                (layer_i>1)?\n");                                                                                                                                                                                                                                                                                               
       fprintf(out1,"                (layer_i==2'b11)?\n");                                                                                                                                                                                                                                                                                        
       fprintf(out1,"                {ctm%d[11:9],vtm%d_%d[11:0],ctm%d[8:6],vtm%d_%d[11:0],ctm%d[5:3],vtm%d_%d[11:0],ctm%d[2:0],vtm%d_%d[11:0]}:\n",i%40,i%40,i/40+4,i%40,i%40,i/40+3,i%40,i%40,i/40+2,i%40,i%40,i/40+1);                                                                                                       
       fprintf(out1,"                {vtm%d_%d[14:12 ],ctm%d[11:9],vtm%d_%d[8:0],vtm%d_%d[14:12 ],ctm%d[8:6],vtm%d_%d[8:0],\n",i%40,i/40+4,i%40,i%40,i/40+4,i%40,i/40+3,i%40,i%40,i/40+3);                                                                                                                                                
       fprintf(out1,"                 vtm%d_%d[14:12 ],ctm%d[5:3],vtm%d_%d[8:0],vtm%d_%d[14:12 ],ctm%d[2:0],vtm%d_%d[8:0]}:\n",i%40,i/40+2,i%40,i%40,i/40+2,i%40,i/40+1,i%40,i%40,i/40+1);                                                                                                                                               
     fprintf(out1,"                (layer_i==2'b01)?\n");                                                                                                                                                                                                                                                                                          
       fprintf(out1,"                {vtm%d_%d[14:9 ],ctm%d[11:9],vtm%d_%d[2:0],vtm%d_%d[14:9 ],ctm%d[8:6],vtm%d_%d[5:0],\n",i%40,i/40+4,i%40,i%40,i/40+4,i%40,i/40+3,i%40,i%40,i/40+3);                                                                                                                                                  
      fprintf(out1, "                 vtm%d_%d[14:9 ],ctm%d[5:3],vtm%d_%d[2:0],vtm%d_%d[14:9 ],ctm%d[2:0],vtm%d_%d[5:0]}:\n",i%40,i/40+2,i%40,i%40,i/40+2,i%40,i/40+1,i%40,i%40,i/40+1);                                                                                                                                                 
     fprintf(out1,"                {vtm%d_%d[14:6],ctm%d[11:9],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[8:6],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[5:3],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm_%d[2:0],vtm%d_%d[2:0]};\n",i%40,i/40+4,i%40,i%40,i/40+4,i%40,i/40+3,i%40,i%40,i/40+3,i%40,i/40+2,i%40,i%40,i/40+2,i%40,i/40+1,i%40,i%40,i/40+1);
*/
 /* 	     fprintf(out1,"wire [119:0] rtm%d=\n",i);                                                                                                                                                                                               
		fprintf(out1,"                (layer_i>1)?\n");                                                                                                                                                                                           
       fprintf(out1,"                (layer_i==2'b11)?\n");                                                                                                                                                                                       
       fprintf(out1,"                {ctm%d[23:21],vtm%d_%d[11:0],ctm%d[20:18],vtm%d_%d[11:0],ctm%d[17:15],vtm%d_%d[11:0],ctm%d[14:12],vtm%d_%d[11:0],\n",i%40,i%40,i/40+8,i%40,i%40,i/40+7,i%40,i%40,i/40+6,i%40,i%40,i/40+5);       
       fprintf(out1,"                 ctm%d[11:9],vtm%d_%d[11:0],ctm%d[8:6],vtm%d_%d[11:0],ctm%d[5:3],vtm%d_%d[11:0],ctm%d[2:0],vtm%d_%d[11:0]}:\n",i%40,i%40,i/40+4,i%40,i%40,i/40+3,i%40,i%40,i/40+2,i%40,i%40,i/40+1);      
       fprintf(out1,"                {vtm%d_%d[14:12 ],ctm%d[23:21],vtm%d_%d[8:0],vtm%d_%d[14:12 ],ctm%d[20:18],vtm%d_%d[8:0],\n",i%40,i/40+8,i%40,i%40,i/40+8,i%40,i/40+7,i%40,i%40,i/40+7);                                             
       fprintf(out1,"                 vtm%d_%d[14:12 ],ctm%d[17:15],vtm%d_%d[8:0],vtm%d_%d[14:12 ],ctm%d[14:12],vtm%d_%d[8:0],\n",i%40,i/40+6,i%40,i%40,i/40+6,i%40,i/40+5,i%40,i%40,i/40+5);                                             
       fprintf(out1,"                 vtm%d_%d[14:12 ],ctm%d[11:9],vtm%d_%d[8:0],vtm%d_%d[14:12 ],ctm%d[8:6],vtm%d_%d[8:0],\n",i%40,i/40+4,i%40,i%40,i/40+4,i%40,i/40+3,i%40,i%40,i/40+3);                                              
       fprintf(out1,"                 vtm%d_%d[14:12 ],ctm%d[5:3],vtm%d_%d[8:0],vtm%d_%d[14:12 ],ctm%d[2:0],vtm%d_%d[8:0]}:\n",i%40,i/40+2,i%40,i%40,i/40+2,i%40,i/40+1,i%40,i%40,i/40+1);                                             
     fprintf(out1,"                (layer_i==2'b01)?\n");                                                                                                                                                                                         
       fprintf(out1,"                {vtm%d_%d[14:9 ],ctm%d[23:21],vtm%d_%d[5:0],vtm%d_%d[14:9 ],ctm%d[20:18],vtm%d_%d[5:0],\n",i%40,i/40+8,i%40,i%40,i/40+8,i%40,i/40+7,i%40,i%40,i/40+7);                                             
       fprintf(out1,"                 vtm%d_%d[14:9 ],ctm%d[17:15],vtm%d_%d[5:0],vtm%d_%d[14:9 ],ctm%d[14:12],vtm%d_%d[5:0],\n",i%40,i/40+6,i%40,i%40,i/40+6,i%40,i/40+5,i%40,i%40,i/40+5);                                             
      fprintf(out1, "                 vtm%d_%d[14:9 ],ctm%d[11:9],vtm%d_%d[5:0],vtm%d_%d[14:9 ],ctm%d[8:6],vtm%d_%d[5:0],\n",i%40,i/40+4,i%40,i%40,i/40+4,i%40,i/40+3,i%40,i%40,i/40+3);                                              
      fprintf(out1, "                 vtm%d_%d[14:9 ],ctm%d[5:3],vtm%d_%d[5:0],vtm%d_%d[14:9 ],ctm%d[2:0],vtm%d_%d[5:0]}:\n",i%40,i/40+2,i%40,i%40,i/40+2,i%40,i/40+1,i%40,i%40,i/40+1);                                             
     fprintf(out1,"                {vtm%d_%d[14:6],ctm%d[23:21],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[20:18],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[17:15],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[14:12],vtm%d_%d[2:0],\n",i%40,i/40+8,i%40,i%40,i/40+8,i%40,i/40+7,i%40,i%40,i/40+7,i%40,i/40+6,i%40,i%40,i/40+6,i%40,i/40+5,i%40,i%40,i/40+5);    
     fprintf(out1,"                 vtm%d_%d[14:6],ctm%d[11:9],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[8:6],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[5:3],vtm%d_%d[2:0],vtm%d_%d[14:6],ctm%d[2:0],vtm%d_%d[2:0]};\n",i%40,i/40+4,i%40,i%40,i/40+4,i%40,i/40+3,i%40,i%40,i/40+3,i%40,i/40+2,i%40,i%40,i/40+2,i%40,i/40+1,i%40,i%40,i/40+1);   
   */
	//  fprintf(out1,"assign mtv%d_%d_w=mtv%d_w[%d:%d];\n",i%40,i/40+1,i%40,3*(i/40)+2,3*(i/40));
	//   if(i%40!=0)fprintf(out1,"wire [14:0] vtm%d_%d=mtv%d_r[%d:%d];\n",i%40,i/40+1,i%40,15*(i/40)+14,15*(i/40)); 
//	   else fprintf(out1,"wire [14:0] vtm%d_%d=mtv%d_r[%d:%d];\n",40,i/40,40,15*(i/40-1)+14,15*(i/40-1)); 
 

	 //   fprintf(out1,"assign summtv%d_%d_w=summtv%d_w[%d:%d];\n",i%40,i/40+1,i%40,4*(i/40)+3,4*(i/40));
/*for(j=0;j<8;j++)
{

fprintf(out1,"assign twosmtv%d[%d:%d]=(mtv%d_w[%d])?{mtv%d_w[%d],~mtv%d_w[%d:%d]+2'b1}:mtv%d_w[%d:%d];\n",i,15*j+5,15*j+3,i,15*j+5,i,15*j+5,i,15*j+4,15*j+3,i,15*j+5,15*j+3);


fprintf(out1,"assign twosmtv%d[%d:%d]=(mtv%d_w[%d])?{mtv%d_w[%d],~mtv%d_w[%d:%d]+2'b1}:mtv%d_w[%d:%d];\n",i,15*j+8,15*j+6,i,15*j+8,i,15*j+8,i,15*j+7,15*j+6,i,15*j+8,15*j+6); 
          

fprintf(out1,"assign twosmtv%d[%d:%d]=(mtv%d_w[%d])?{mtv%d_w[%d],~mtv%d_w[%d:%d]+2'b1}:mtv%d_w[%d:%d];\n",i,15*j+11,15*j+9,i,15*j+11,i,15*j+11,i,15*j+10,15*j+9,i,15*j+11,15*j+9);
         

fprintf(out1,"assign twosmtv%d[%d:%d]=(mtv%d_w[%d])?{mtv%d_w[%d],~mtv%d_w[%d:%d]+2'b1}:mtv%d_w[%d:%d];\n",i,15*j+14,15*j+12,i,15*j+14,i,15*j+14,i,15*j+13,15*j+12,i,15*j+14,15*j+12);
         

fprintf(out1,"assign twosmtv%d[%d:%d]=(mtv%d_w[%d])?{mtv%d_w[%d],~mtv%d_w[%d:%d]+2'b1}:mtv%d_w[%d:%d];\n",i,15*j+2,15*j,i,15*j+2,i,15*j+2,i,15*j+1,15*j,i,15*j+2,15*j);
           



}
*/
	   
	   for(j=0;j<8;j++)
	   {
      //  fprintf(out1,"sum%d_w[%d:%d] = mtv%d_w[%d:%d] + mtv%d_w[%d:%d] + mtv%d_w[%d:%d] + mtv%d_w[%d:%d] + mtv%d_w[%d:%d];\n",i,4*j+3,4*j,i,15*j+14,15*j+12,i,15*j+11,15*j+9,i,15*j+8,15*j+6,i,15*j+5,15*j+3,i,15*j+2,15*j);
	   // fprintf(out1,"stv%d_w[%d:%d]=(layer_i<2)?(layer_i==2'b00)?mtv%d_w[%d:%d]:mtv%d_w[%d:%d]:(layer_i==2'b10)?mtv%d_w[%d:%d]:mtv%d_w[%d:%d];\n",i,3*j+2,3*j,i,15*j+5,15*j+3,i,15*j+8,15*j+6,i,15*j+11,15*j+9,i,15*j+14,15*j+12);  

	   }
/*	fprintf(out1,"memblock mem%d_b(\n",i);
	fprintf(out1,"           .address(addr%db_w),\n",i);
	fprintf(out1,"         .clken(mem_en),\n");
	fprintf(out1,"         .clock(CLK),\n");
	fprintf(out1,"         .data(rtm%db_w),\n",i);
	fprintf(out1,"         .wren(memwrb%d_en),\n",i);
	fprintf(out1,"         .q(mtv%db)\n",i);
	fprintf(out1,"         );\n");*/
/*	   fprintf(out1,"memblock32x120 mem%d_a(\n",i);
	fprintf(out1,"           .addra(addr%d_w),\n",i);
	fprintf(out1,"         .ena(dec_en),\n");
	fprintf(out1,"         .clka(CLK),\n");
	fprintf(out1,"         .dina(rtm%da_w),\n",i);
	fprintf(out1,"         .wea(memwra%d_en),\n",i);
	fprintf(out1,"         .douta(mtv%da)\n",i);
	fprintf(out1,"         );\n");*/
// if(i%40!=0) fprintf(out1,"wire [2:0] chdata%d_%d_w =chdata%d_w[%d:%d];\n",i%40,i/40+1,i/40+1,3*(i/40)+2,3*(i/40));
// else  fprintf(out1,"wire [2:0] chdata%d_%d_w =chdata%d_w[%d:%d];\n",40,i/40,i/40,3*(i/40-1)+2,3*(i/40-1));
/*fprintf(out1,"chdata%d[2:0]= mtv%d_w[2:0];\n",i,i);
fprintf(out1,"chdata%d[5:3]= mtv%d_w[17:15];\n",i,i); 
fprintf(out1,"chdata%d[8:6]= mtv%d_w[32:30];\n",i,i); 
fprintf(out1,"chdata%d[11:9]=mtv%d_w[47:45];\n",i,i); */
/*
fprintf(out1,"chdata%d[14:12]= mtv%d_w[62:60];\n",i,i);
fprintf(out1,"chdata%d[17:15]= mtv%d_w[77:75];\n",i,i); 
fprintf(out1,"chdata%d[20:18]= mtv%d_w[92:90];\n",i,i); 
fprintf(out1,"chdata%d[23:21]=mtv%d_w[107:105];\n",i,i);  */
/*
fprintf(out1,"      6'd%d:\n",i+32);
fprintf(out1,"      begin\n");
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[2:0];\n"  ,8*(i%4)+0,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[5:3];\n"  ,8*(i%4)+1,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[8:6];\n"  ,8*(i%4)+2,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[11:9];\n" ,8*(i%4)+3,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[14:12];\n",8*(i%4)+4,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[17:15];\n",8*(i%4)+5,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[20:18];\n",8*(i%4)+6,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r2 [%d][%d:%d]<=chdata[23:21];\n",8*(i%4)+7,(i/4)*3+2,(i/4)*3);                 
fprintf(out1,"      end\n");
*/
/*
fprintf(out1,"      6'd%d:\n",i);
fprintf(out1,"      begin\n");
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[2:0];\n"  ,8*(i%4)+0,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[5:3];\n"  ,8*(i%4)+1,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[8:6];\n"  ,8*(i%4)+2,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[11:9];\n" ,8*(i%4)+3,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[14:12];\n",8*(i%4)+4,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[17:15];\n",8*(i%4)+5,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[20:18];\n",8*(i%4)+6,(i/4)*3+2,(i/4)*3);
fprintf(out1,"        chdata_r1 [%d][%d:%d]<=chdata[23:21];\n",8*(i%4)+7,(i/4)*3+2,(i/4)*3);                 
fprintf(out1,"      end\n");
*/
   }
    fclose(out1);
	system("pause");
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
	  BitNodeIdx = QC436->CN[i].BitNodeIdx[j] ;
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
  double Mji;
  double TotalPhi;
  int TotalSign;
//**Compute each circulant block of check node than go to bit node process promptly 
  for(i=0 ; i<CheckNodeNum ; i++)/***************/
  {
	TotalPhi = 0 ;
	TotalSign = 1 ;
	
	//Compute total phi and sign for ALL bit node
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx = QC436->CN[i].BitNodeIdx[j] ;
	  ExtArrayIdx = QC436->BN[BitNodeIdx].index ;
	  Mji = QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx] ;
	  TotalPhi += PhiFunction(Mji) ;
	  TotalSign *= SignFunction(Mji) ;
    }

    // exclude jth branch and compute the extrinsic information
	for(j=0 ; j<CheckNodeDeg ; j++)
	{
	  BitNodeIdx = QC436->CN[i].BitNodeIdx[j] ;
	  ExtArrayIdx = QC436->BN[BitNodeIdx].index ;
	  QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx] =
		  (TotalSign * SignFunction(QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx]))  // total sign exclude self sign
		  * PhiFunction(TotalPhi - PhiFunction(QC436->BN[BitNodeIdx].ExtInf[ExtArrayIdx])) ; // total phi exclude self phi
	  QC436->BN[BitNodeIdx].index = (QC436->BN[BitNodeIdx].index+1) % BitNodeDeg ; // index for extrinsic array +1
	}
  }
}

void BitNodeProcessing(Ldpc *QC436)
{
  int i,j ;
  double TotalLlr ;
  // Update the bit node related to the check node process, however some bit node is unnecessary to update.
  for(i=0 ; i<BitNodeNum ; i++)
  {
	TotalLlr = QC436->BN[i].RecLlr ;
	for(j=0 ; j<BitNodeDeg ; j++) 
	  TotalLlr += QC436->BN[i].ExtInf[j];
    
	QC436->BN[i].TempHD = (TotalLlr < 0) ? 1 : 0 ;

	for(j=0 ; j<BitNodeDeg ; j++)
	  QC436->BN[i].ExtInf[j] = TotalLlr - QC436->BN[i].ExtInf[j];
  }
}

Ldpc QC436;
void main(void)
{
  int i,j;
  long seed = -1;
  FILE *out1;
  double SNR;         /* in dB*/
  double var;
  int IterNum;
  int BlockNum;
  int ErrorFlag, ErrorBlockNumber, ErrorBits, RawErrorBits;
  int ParityCheckResult;

  FFGenHMatrix(&QC436) ;


  out1=fopen("FinalResult.txt","a+t");
  fprintf(out1, "\nQC-LDPC with %d Bit Nodes and %d Check Nodes and Circulant Size = %d", BitNodeNum, CheckNodeNum, CirculantSize);
  fprintf(out1, "\nSum-Product Decoding Algorithm with double input and Maximum %d Iterations", MaxIterNum);
  fprintf(out1, "\nBit Node Degree = %d, Check Node Degree = %d, with Random H Matrix with Seed = %d", BitNodeDeg, CheckNodeDeg, InitialMatrixSeed);
  fprintf(out1, "\n\n") ;
  fclose(out1);

  //Assume All Zero Codeword Is Transmitted so that we can ignore LDPC encoder at this stage
  for(SNR=3.0 ; SNR<=6.0 ; SNR+=0.1)
  {
	printf("\nSNR=%f", SNR) ;
    var=(double)((BitNodeNum)/(BitNodeNum-CheckNodeNum)/2.0*pow(10,-(double)SNR/10));

    ErrorBlockNumber = 0 ;
	ErrorBits = 0 ;
	RawErrorBits = 0 ;
	for(BlockNum=0 ; ErrorBlockNumber<20 ; BlockNum++)
	//for(BlockNum=0 ; BlockNum<1000 ; BlockNum++)
	{
	  if(BlockNum % 100 == 0)
	  {
	    printf(".") ;
		out1=fopen("Result.txt","w");
        fprintf(out1, "\nSNR:%f, RawBER=%d/%d=%e, BER=%d/%d=%e, BLER=%d/%d=%e\n" , SNR, RawErrorBits, BitNodeNum*BlockNum, (double)RawErrorBits/(BitNodeNum*BlockNum), 
		  ErrorBits, BitNodeNum*BlockNum, (double)ErrorBits/(BitNodeNum*BlockNum), ErrorBlockNumber, BlockNum ,(double)ErrorBlockNumber/(BlockNum) );
        fclose(out1);
	  }
	  
	  // Channel Noise and Initialization
	  for(i=0 ; i<BitNodeNum ; i++)
	  {
        QC436.BN[i].RecLlr = (double)(2 * (1 + sqrt( var)*gaussian_noise(&seed)) / var) ;
		QC436.BN[i].TempHD = (QC436.BN[i].RecLlr <0 ) ? 1 : 0 ;
	    for(j=0 ; j<BitNodeDeg ; j++)
	      QC436.BN[i].ExtInf[j] = QC436.BN[i].RecLlr ;
		QC436.BN[i].index = 0 ;     
		
		if(QC436.BN[i].TempHD == 1)
	      RawErrorBits ++ ;
	  }

	  IterNum = 0 ;
	  ParityCheckResult = ParityCheck(&QC436) ;
	  for(IterNum=0 ; IterNum<MaxIterNum && ParityCheckResult==1 ; IterNum++)
	  {
        
        CheckNodeProcessing(i,&QC436) ;
        BitNodeProcessing(&QC436) ;
		
		ParityCheckResult = ParityCheck(&QC436) ;
	  } // End of for(IterNum=0 ; IterNum<MaxIterNum ; IterNum++)

	  ErrorFlag=0;
      for(i=0 ; i<BitNodeNum ; i++)
	  {
        if(QC436.BN[i].TempHD != 0)
		{
	      ErrorBits++;
          ErrorFlag=1;
		}
	  }
      if(ErrorFlag==1)
	    ErrorBlockNumber++;
 
    } // End of for(BlockNum=0 ; BlockNum<100 ; BlockNum++)

    out1=fopen("FinalResult.txt","a+t");
    fprintf(out1, "\nSNR:%f, RawBER=%d/%d=%e, BER=%d/%d=%e, BLER=%d/%d=%e\n" , SNR, RawErrorBits, BitNodeNum*BlockNum, (double)RawErrorBits/(BitNodeNum*BlockNum), 
		ErrorBits, BitNodeNum*BlockNum, (double)ErrorBits/((double)BitNodeNum*BlockNum), ErrorBlockNumber, BlockNum ,(double)ErrorBlockNumber/(BlockNum) );
    fclose(out1);

  } // End of for(SNR=1.0 ; SNR<=5.0 ; SNR+=1.0)

}