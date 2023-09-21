//----------------------------------------------------------------------------------------
//
// $Author: haschou $
// $Date: 2012/02/06 $
// $Revision: 0.1 $
//----------------------------------------------------------------------------------------
// Project: LDPC Decoder
//            Testbench of LDPC decoder 48bit input/ 40*8bit output
//            Input code data from file codeword.txt
//            Test algorithm behavior for each single block
//----------------------------------------------------------------------------------------
//sqrt(Var)/Var
//0111111111
//1000000101
//1000001011
//1000010001
//1000010111
//1000011101
//1000100011
//1000101010
//1000110000
//1000110111
`timescale 1ns/10ps
`define CLKPERIOD 16
module LDPCDECtb;
reg       CLK,RESET_N,start_ldpc,start_Ldata,start_encode;
reg [15:0] data_b;
reg [47:0] data_in;
reg [47:0] data_in_r;
wire [47:0] data_in_w;
reg [15:0] data_out;
reg [47:0] data;
reg [9:0] gaussnoise;
reg [9:0] gnacc;
reg [19:0] gnoise;
reg [3:0]  data_gn;
wire [3:0]  hb_count;
wire [15:0] hardbit,parity_r;
wire enga_en,hb_vlid,encode_fin;
wire [12:0] ch_counta,ch_countb;
wire [4:0] maxiter=5'd8;
integer i,j,k,n,m,fptr_o,berr,inx,e,p,l,err;
parameter block_size=10240;
parameter num_block=9;
parameter data_leng=num_block*block_size;
parameter idle_time=256+16;
parameter idle_leng=(num_block-1)*idle_time;
parameter iteration_time=768*8*16;
parameter hardbitout_time=320;
reg  [14:0] err_loc [0:1000];
reg  [47:0] packetmem [0:num_block*block_size];
reg  [47:0] packetmem1 [0:num_block*block_size];
reg  [9:0] packetmemgn [0:1023];
reg  [14:0] err_loc1;
LDPCg5p16codec LDPCg5p16codec(
                              .CLK(CLK),
                              .RESET_N(RESET_N),
                              .start_encode(start_encode),
                              .maxiter(maxiter),
                              .data_in_tx(data_b),
                              .data_in(data_in),
                              .encode_fin(encode_fin),
                              .parity_r(parity_r),
                              .start_ldpc(start_ldpc),
                              .start_Ldata(start_Ldata),
                              .hb_count(hb_count),
                              .hb_vlid(hb_vlid),
                              .hardbit(hardbit),
                              .enga_en(enga_en),
                              .dec_fin(dec_fin),
                              .dec_fail(dec_fail),
                              .ch_counta(ch_counta),.ch_countb(ch_countb));

initial
begin
 CLK=1'b1;
 RESET_N=1'b0;
 start_ldpc=0;
 start_Ldata=0;
 start_encode=0;
 m=0;
 i=0;
 berr=0;
 k=0;
 j=0;
 e=0;
 p=0;
 l=0;
 err=0;
end

initial forever #(`CLKPERIOD) CLK = ~CLK;

initial
begin

$fsdbDumpfile("LDPCcodec_rtl.fsdb");
$fsdbDumpvars;
$readmemb("codeword.txt",packetmem);
$readmemb("codeword1.txt",packetmem1); //correct codeword

$readmemb("err_loc.txt",err_loc);
$readmemb("GaussNoise.txt",packetmemgn);
 fptr_o = $fopen("hwout.dat");
  if(!fptr_o) begin
    $display("Can not write file!");
    $finish;
  end

 #160 RESET_N=1'b0;
 RESET_N=1'b1;


 start_ldpc=1;
  start_Ldata=1;
      @(posedge CLK);
  for (j=0;j<=block_size+(num_block+1)*iteration_time;j=j+16)
  begin
          // Input data
    if (j<=block_size)
    begin

      start_encode=1;
      data_in_r=packetmem[i];
      data_in={data_in_r[2:0],data_in_r[5:3],data_in_r[8:6],data_in_r[11:9],
               data_in_r[14:12],data_in_r[17:15],data_in_r[20:18],data_in_r[23:21],
               data_in_r[26:24],data_in_r[29:27],data_in_r[32:30],data_in_r[35:33],
               data_in_r[38:36],data_in_r[41:39],data_in_r[44:42],data_in_r[47:45]};
      i=i+1;
      @(posedge CLK);
    end
    else if(j>block_size && j<=block_size+idle_time)
    begin
      start_Ldata=1;
      start_encode=0;

      @(posedge CLK);
    end
    else if(j>block_size+idle_time && j<=2*block_size+idle_time)
    begin
      e=1;
      start_Ldata=1;
      start_encode=1;
      data_in_r=packetmem[i];
      data_in={data_in_r[2:0],data_in_r[5:3],data_in_r[8:6],data_in_r[11:9],
               data_in_r[14:12],data_in_r[17:15],data_in_r[20:18],data_in_r[23:21],
               data_in_r[26:24],data_in_r[29:27],data_in_r[32:30],data_in_r[35:33],
               data_in_r[38:36],data_in_r[41:39],data_in_r[44:42],data_in_r[47:45]};
      i=i+1;
       @(posedge CLK);
    end
    else if(j>2*block_size+idle_time && j<=2*block_size+2*idle_time)
    begin
      start_Ldata=1;
      start_encode=0;

      @(posedge CLK);
    end
    else if(j==block_size+idle_time+(m+1)*iteration_time)
    begin
      $display("%dth block %dbit error %d errloc mis",m,berr,err);
           if(e<2)e=e+1;
           else e=0;
           m=m+1;
        @(posedge CLK);
    end
    else if(j>block_size + idle_time + m*iteration_time && j<=2*block_size + idle_time + m*iteration_time)
    begin
      start_Ldata=1;
      start_encode=1;
      data_in_r=packetmem[i];
      data_in={data_in_r[2:0],data_in_r[5:3],data_in_r[8:6],data_in_r[11:9],
               data_in_r[14:12],data_in_r[17:15],data_in_r[20:18],data_in_r[23:21],
               data_in_r[26:24],data_in_r[29:27],data_in_r[32:30],data_in_r[35:33],
               data_in_r[38:36],data_in_r[41:39],data_in_r[44:42],data_in_r[47:45]};
      i=i+1;
      @(posedge CLK);
    end
    else if(j>2*block_size+idle_time+ m*iteration_time && j<=2*block_size+2*idle_time+ m*iteration_time)
    begin
      start_Ldata=1;
      start_encode=0;

      @(posedge CLK);
    end
    else
    begin
      start_Ldata=0;
       start_encode=0;
      @(posedge CLK);
    //  i=0;
    end

         //     Recieve data
    if(hb_vlid)
    begin
      data=packetmem1[k/16];
       for (n=0;n<16;n=n+1)
       begin
         $fwrite(fptr_o,"%d\n",hardbit[n]);
         if( hardbit[n] != data[47-n*3])
         begin

            err_loc1=err_loc[berr];
             berr=berr+1;
            if(err_loc1<k-(m-1)*10240||err_loc1>k+16-(m-1)*10240)
            err=err+1;
         end

       end
       k=k+16;
    end
    else
    begin
 //     k=0;
      p=0;
    end
  end
   $display("%d block Input Fin \n",m);

  $finish;
end
endmodule
