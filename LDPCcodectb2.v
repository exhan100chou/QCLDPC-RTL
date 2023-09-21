//----------------------------------------------------------------------------------------
//
// $Author: haschou $
// $Date: 2012/02/06 $
// $Revision: 0.1 $
//----------------------------------------------------------------------------------------
// Project: LDPC Decoder
//            Testbench of LDPC decoder 48bit input/ 40*8bit output
//            Gaussian noise gen can run desired block number
//            Test the BER performance
//----------------------------------------------------------------------------------------

`timescale 1ns/10ps
`define CLKPERIOD 16
module LDPCDECtb;
reg       CLK,RESET_N,start_ldpc,start_Ldata,start_encode;
reg [15:0] data_b;
wire [4:0] maxiter=5'd8;
wire [47:0] data_in;
wire [47:0] data_in_w;
wire [63:0]  data_in_gn;
reg [15:0] data,data_out;
reg signed [9:0] gaussnoise;
reg [9:0] gnacc;
reg signed [19:0] gnoise;
reg [47:0]  data_gn;
wire [3:0]  hb_count;
wire [15:0] hardbit,parity_r;
wire enga_en,hb_vlid,encode_fin;
wire [12:0] ch_counta,ch_countb;
wire [9:0] var1 =10'b01011111110;
wire [9:0] var2 =10'b01011111101;
wire [9:0] var3 =10'b01011111100;
wire [9:0] var4 =10'b01011111011;
wire [9:0] var5 =10'b01011111010;
wire [9:0] var6 =10'b01011111000;
wire [9:0] var7 =10'b01011110111;
wire [9:0] var8 =10'b01011110110;
wire [9:0] var9 =10'b01011110101;
wire [9:0] var10=10'b01011110011;
reg signed [9:0] var;
integer i,j,k,n,m,fptr_o,berr,inx,e,p,l,varinx,RBER;
real BER,A,B,rbber,C;
parameter block_size=9216;
parameter num_block=3;
parameter data_leng=num_block*block_size;
parameter idle_time=1024+256+16;
parameter idle_leng=(num_block-1)*idle_time;
parameter iteration_time=768*8*16;
parameter hardbitout_time=320;
reg  [15:0] packetmem1 [0:640];
reg  [15:0] packetmem2 [0:640];
reg  [15:0] packetmem3 [0:640];
reg  [9:0] packetmemgn [0:1023];
wire [3:0] addr_count;
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
                              .addr_count(addr_count),
                              .hardbit(hardbit),
                              .enga_en(enga_en),
                              .dec_fin(dec_fin),
                              .dec_fail(dec_fail),
                              .ch_counta(ch_counta),.ch_countb(ch_countb));                                                                    
assign data_in_w[47:45]=(start_encode)? data_b[ 0] ? {data_b[ 0],2'b01}:{data_b[ 0],2'b11}:(encode_fin)? parity_r[ 0]? {parity_r[ 0],2'b01}:{parity_r[ 0],2'b11}: 3'd0;
assign data_in_w[44:42]=(start_encode)? data_b[ 1] ? {data_b[ 1],2'b01}:{data_b[ 1],2'b11}:(encode_fin)? parity_r[ 1]? {parity_r[ 1],2'b01}:{parity_r[ 1],2'b11}: 3'd0;
assign data_in_w[41:39]=(start_encode)? data_b[ 2] ? {data_b[ 2],2'b01}:{data_b[ 2],2'b11}:(encode_fin)? parity_r[ 2]? {parity_r[ 2],2'b01}:{parity_r[ 2],2'b11}: 3'd0;
assign data_in_w[38:36]=(start_encode)? data_b[ 3] ? {data_b[ 3],2'b01}:{data_b[ 3],2'b11}:(encode_fin)? parity_r[ 3]? {parity_r[ 3],2'b01}:{parity_r[ 3],2'b11}: 3'd0;
assign data_in_w[35:33]=(start_encode)? data_b[ 4] ? {data_b[ 4],2'b01}:{data_b[ 4],2'b11}:(encode_fin)? parity_r[ 4]? {parity_r[ 4],2'b01}:{parity_r[ 4],2'b11}: 3'd0;
assign data_in_w[32:30]=(start_encode)? data_b[ 5] ? {data_b[ 5],2'b01}:{data_b[ 5],2'b11}:(encode_fin)? parity_r[ 5]? {parity_r[ 5],2'b01}:{parity_r[ 5],2'b11}: 3'd0;
assign data_in_w[29:27]=(start_encode)? data_b[ 6] ? {data_b[ 6],2'b01}:{data_b[ 6],2'b11}:(encode_fin)? parity_r[ 6]? {parity_r[ 6],2'b01}:{parity_r[ 6],2'b11}: 3'd0;
assign data_in_w[26:24]=(start_encode)? data_b[ 7] ? {data_b[ 7],2'b01}:{data_b[ 7],2'b11}:(encode_fin)? parity_r[ 7]? {parity_r[ 7],2'b01}:{parity_r[ 7],2'b11}: 3'd0;
assign data_in_w[23:21]=(start_encode)? data_b[ 8] ? {data_b[ 8],2'b01}:{data_b[ 8],2'b11}:(encode_fin)? parity_r[ 8]? {parity_r[ 8],2'b01}:{parity_r[ 8],2'b11}: 3'd0;
assign data_in_w[20:18]=(start_encode)? data_b[ 9] ? {data_b[ 9],2'b01}:{data_b[ 9],2'b11}:(encode_fin)? parity_r[ 9]? {parity_r[ 9],2'b01}:{parity_r[ 9],2'b11}: 3'd0;
assign data_in_w[17:15]=(start_encode)? data_b[10] ? {data_b[10],2'b01}:{data_b[10],2'b11}:(encode_fin)? parity_r[10]? {parity_r[10],2'b01}:{parity_r[10],2'b11}: 3'd0;
assign data_in_w[14:12]=(start_encode)? data_b[11] ? {data_b[11],2'b01}:{data_b[11],2'b11}:(encode_fin)? parity_r[11]? {parity_r[11],2'b01}:{parity_r[11],2'b11}: 3'd0;
assign data_in_w[11: 9]=(start_encode)? data_b[12] ? {data_b[12],2'b01}:{data_b[12],2'b11}:(encode_fin)? parity_r[12]? {parity_r[12],2'b01}:{parity_r[12],2'b11}: 3'd0;
assign data_in_w[ 8: 6]=(start_encode)? data_b[13] ? {data_b[13],2'b01}:{data_b[13],2'b11}:(encode_fin)? parity_r[13]? {parity_r[13],2'b01}:{parity_r[13],2'b11}: 3'd0;
assign data_in_w[ 5: 3]=(start_encode)? data_b[14] ? {data_b[14],2'b01}:{data_b[14],2'b11}:(encode_fin)? parity_r[14]? {parity_r[14],2'b01}:{parity_r[14],2'b11}: 3'd0;
assign data_in_w[ 2: 0]=(start_encode)? data_b[15] ? {data_b[15],2'b01}:{data_b[15],2'b11}:(encode_fin)? parity_r[15]? {parity_r[15],2'b01}:{parity_r[15],2'b11}: 3'd0;

wire signed [2:0]
data_in1_w =data_in_w[47:45],
data_in2_w =data_in_w[44:42],
data_in3_w =data_in_w[41:39],
data_in4_w =data_in_w[38:36],
data_in5_w =data_in_w[35:33],
data_in6_w =data_in_w[32:30],
data_in7_w =data_in_w[29:27],
data_in8_w =data_in_w[26:24],
data_in9_w =data_in_w[23:21],
data_in10_w=data_in_w[20:18],
data_in11_w=data_in_w[17:15],
data_in12_w=data_in_w[14:12],
data_in13_w=data_in_w[11: 9],
data_in14_w=data_in_w[ 8: 6],
data_in15_w=data_in_w[ 5: 3],
data_in16_w=data_in_w[ 2: 0],
data_gn1 =data_gn[47:45],
data_gn2 =data_gn[44:42],
data_gn3 =data_gn[41:39],
data_gn4 =data_gn[38:36],
data_gn5 =data_gn[35:33],
data_gn6 =data_gn[32:30],
data_gn7 =data_gn[29:27],
data_gn8 =data_gn[26:24],
data_gn9 =data_gn[23:21],
data_gn10=data_gn[20:18],
data_gn11=data_gn[17:15],
data_gn12=data_gn[14:12],
data_gn13=data_gn[11: 9],
data_gn14=data_gn[ 8: 6],
data_gn15=data_gn[ 5: 3],
data_gn16=data_gn[ 2: 0];


wire signed [3:0]
 data_in_gn1_w ={data_in1_w [2],data_in1_w }+{data_gn1 [2],data_gn1 },
 data_in_gn2_w ={data_in2_w [2],data_in2_w }+{data_gn2 [2],data_gn2 },
 data_in_gn3_w ={data_in3_w [2],data_in3_w }+{data_gn3 [2],data_gn3 },
 data_in_gn4_w ={data_in4_w [2],data_in4_w }+{data_gn4 [2],data_gn4 },
 data_in_gn5_w ={data_in5_w [2],data_in5_w }+{data_gn5 [2],data_gn5 },
 data_in_gn6_w ={data_in6_w [2],data_in6_w }+{data_gn6 [2],data_gn6 },
 data_in_gn7_w ={data_in7_w [2],data_in7_w }+{data_gn7 [2],data_gn7 },
 data_in_gn8_w ={data_in8_w [2],data_in8_w }+{data_gn8 [2],data_gn8 },
 data_in_gn9_w ={data_in9_w [2],data_in9_w }+{data_gn9 [2],data_gn9 },
 data_in_gn10_w={data_in10_w[2],data_in10_w}+{data_gn10[2],data_gn10},
 data_in_gn11_w={data_in11_w[2],data_in11_w}+{data_gn11[2],data_gn11},
 data_in_gn12_w={data_in12_w[2],data_in12_w}+{data_gn12[2],data_gn12},
 data_in_gn13_w={data_in13_w[2],data_in13_w}+{data_gn13[2],data_gn13},
 data_in_gn14_w={data_in14_w[2],data_in14_w}+{data_gn14[2],data_gn14},
 data_in_gn15_w={data_in15_w[2],data_in15_w}+{data_gn15[2],data_gn15},
 data_in_gn16_w={data_in16_w[2],data_in16_w}+{data_gn16[2],data_gn16};



assign data_in[47:45]=(data_in_gn1_w [2] ^ data_in_gn1_w [3])? {data_in_gn1_w [3],2'b11} : {data_in_gn1_w [3],data_in_gn1_w [1:0]};
assign data_in[44:42]=(data_in_gn2_w [2] ^ data_in_gn2_w [3])? {data_in_gn2_w [3],2'b11} : {data_in_gn2_w [3],data_in_gn2_w [1:0]};
assign data_in[41:39]=(data_in_gn3_w [2] ^ data_in_gn3_w [3])? {data_in_gn3_w [3],2'b11} : {data_in_gn3_w [3],data_in_gn3_w [1:0]};
assign data_in[38:36]=(data_in_gn4_w [2] ^ data_in_gn4_w [3])? {data_in_gn4_w [3],2'b11} : {data_in_gn4_w [3],data_in_gn4_w [1:0]};
assign data_in[35:33]=(data_in_gn5_w [2] ^ data_in_gn5_w [3])? {data_in_gn5_w [3],2'b11} : {data_in_gn5_w [3],data_in_gn5_w [1:0]};
assign data_in[32:30]=(data_in_gn6_w [2] ^ data_in_gn6_w [3])? {data_in_gn6_w [3],2'b11} : {data_in_gn6_w [3],data_in_gn6_w [1:0]};
assign data_in[29:27]=(data_in_gn7_w [2] ^ data_in_gn7_w [3])? {data_in_gn7_w [3],2'b11} : {data_in_gn7_w [3],data_in_gn7_w [1:0]};
assign data_in[26:24]=(data_in_gn8_w [2] ^ data_in_gn8_w [3])? {data_in_gn8_w [3],2'b11} : {data_in_gn8_w [3],data_in_gn8_w [1:0]};
assign data_in[23:21]=(data_in_gn9_w [2] ^ data_in_gn9_w [3])? {data_in_gn9_w [3],2'b11} : {data_in_gn9_w [3],data_in_gn9_w [1:0]};
assign data_in[20:18]=(data_in_gn10_w[2] ^ data_in_gn10_w[3])? {data_in_gn10_w[3],2'b11} : {data_in_gn10_w[3],data_in_gn10_w[1:0]};
assign data_in[17:15]=(data_in_gn11_w[2] ^ data_in_gn11_w[3])? {data_in_gn11_w[3],2'b11} : {data_in_gn11_w[3],data_in_gn11_w[1:0]};
assign data_in[14:12]=(data_in_gn12_w[2] ^ data_in_gn12_w[3])? {data_in_gn12_w[3],2'b11} : {data_in_gn12_w[3],data_in_gn12_w[1:0]};
assign data_in[11: 9]=(data_in_gn13_w[2] ^ data_in_gn13_w[3])? {data_in_gn13_w[3],2'b11} : {data_in_gn13_w[3],data_in_gn13_w[1:0]};
assign data_in[ 8: 6]=(data_in_gn14_w[2] ^ data_in_gn14_w[3])? {data_in_gn14_w[3],2'b11} : {data_in_gn14_w[3],data_in_gn14_w[1:0]};
assign data_in[ 5: 3]=(data_in_gn15_w[2] ^ data_in_gn15_w[3])? {data_in_gn15_w[3],2'b11} : {data_in_gn15_w[3],data_in_gn15_w[1:0]};
assign data_in[ 2: 0]=(data_in_gn16_w[2] ^ data_in_gn16_w[3])? {data_in_gn16_w[3],2'b11} : {data_in_gn16_w[3],data_in_gn16_w[1:0]};

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
 RBER=0;
end

initial forever #(`CLKPERIOD) CLK = ~CLK;

initial
begin
if(num_block<5)
begin
  $fsdbDumpfile("LDPCcodec_rtl1.fsdb");
  $fsdbDumpvars;
end
//$readmemb("message.txt",packetmem1);
$readmemb("GaussNoise.txt",packetmemgn);
 fptr_o = $fopen("hwout.dat");
  if(!fptr_o) begin
    $display("Can not write file!");
    $finish;
  end


  for (varinx=1; varinx<10 ; varinx=varinx+1)
  begin
 RESET_N=1'b0;
 #160 RESET_N=1'b0;
 RESET_N=1'b1;
 start_ldpc=1;
      e=0;
     case(varinx)
     0:var=var1;
     1:var=var2;
     2:var=var3;
     3:var=var4;
     4:var=var5;
     5:var=var6;
     6:var=var7;
     7:var=var8;
     8:var=var9;
     9:var=var10;
     endcase
  for (j=0;j<=block_size+(num_block+1)*iteration_time;j=j+16)
    begin
          // Input data
    if (j<=block_size)
    begin
      start_Ldata=1;
      start_encode=1;
               @(posedge CLK);
     //  data_b=packetmem1[i];
         data_b=$random%32768;
        for(l=0;l<16;l=l+1)
       begin
        gnacc=$random%32768;
        gaussnoise=packetmemgn[gnacc];
        gnoise=gaussnoise*var;
        case(l)
          0 :data_gn[47:45]={gnoise[19],gnoise[16:15]};
          1 :data_gn[44:42]={gnoise[19],gnoise[16:15]};
          2 :data_gn[41:39]={gnoise[19],gnoise[16:15]};
          3 :data_gn[38:36]={gnoise[19],gnoise[16:15]};
          4 :data_gn[35:33]={gnoise[19],gnoise[16:15]};
          5 :data_gn[32:30]={gnoise[19],gnoise[16:15]};
          6 :data_gn[29:27]={gnoise[19],gnoise[16:15]};
          7 :data_gn[26:24]={gnoise[19],gnoise[16:15]};
          8 :data_gn[23:21]={gnoise[19],gnoise[16:15]};
          9 :data_gn[20:18]={gnoise[19],gnoise[16:15]};
          10:data_gn[17:15]={gnoise[19],gnoise[16:15]};
          11:data_gn[14:12]={gnoise[19],gnoise[16:15]};
          12:data_gn[11: 9]={gnoise[19],gnoise[16:15]};
          13:data_gn[ 8: 6]={gnoise[19],gnoise[16:15]};
          14:data_gn[ 5: 3]={gnoise[19],gnoise[16:15]};
          15:data_gn[ 2: 0]={gnoise[19],gnoise[16:15]};
        endcase
         case(l)
           0 :  if( data_in1_w [2] ^ data_in_gn1_w [3] )RBER=RBER+1;
           1 :  if( data_in2_w [2] ^ data_in_gn2_w [3] )RBER=RBER+1;
           2 :  if( data_in3_w [2] ^ data_in_gn3_w [3] )RBER=RBER+1;
           3 :  if( data_in4_w [2] ^ data_in_gn4_w [3] )RBER=RBER+1;
           4 :  if( data_in5_w [2] ^ data_in_gn5_w [3] )RBER=RBER+1;
           5 :  if( data_in6_w [2] ^ data_in_gn6_w [3] )RBER=RBER+1;
           6 :  if( data_in7_w [2] ^ data_in_gn7_w [3] )RBER=RBER+1;
           7 :  if( data_in8_w [2] ^ data_in_gn8_w [3] )RBER=RBER+1;
           8 :  if( data_in9_w [2] ^ data_in_gn9_w [3] )RBER=RBER+1;
           9 :  if( data_in10_w[2] ^ data_in_gn10_w[3] )RBER=RBER+1;
           10:  if( data_in11_w[2] ^ data_in_gn11_w[3] )RBER=RBER+1;
           11:  if( data_in12_w[2] ^ data_in_gn12_w[3] )RBER=RBER+1;
           12:  if( data_in13_w[2] ^ data_in_gn13_w[3] )RBER=RBER+1;
           13:  if( data_in14_w[2] ^ data_in_gn14_w[3] )RBER=RBER+1;
           14:  if( data_in15_w[2] ^ data_in_gn15_w[3] )RBER=RBER+1;
           15:  if( data_in16_w[2] ^ data_in_gn16_w[3] )RBER=RBER+1;
        endcase
       end
    packetmem1[i]=data_b;
      i=i+1;
    end
    else if(j>block_size && j<=block_size+idle_time)
    begin
        for(l=0;l<16;l=l+1)
       begin
        gnacc=$random%1024;
        gaussnoise=packetmemgn[gnacc];
        gnoise=gaussnoise*var;
        case(l)
          0 :data_gn[47:45]={gnoise[19],gnoise[16:15]};
          1 :data_gn[44:42]={gnoise[19],gnoise[16:15]};
          2 :data_gn[41:39]={gnoise[19],gnoise[16:15]};
          3 :data_gn[38:36]={gnoise[19],gnoise[16:15]};
          4 :data_gn[35:33]={gnoise[19],gnoise[16:15]};
          5 :data_gn[32:30]={gnoise[19],gnoise[16:15]};
          6 :data_gn[29:27]={gnoise[19],gnoise[16:15]};
          7 :data_gn[26:24]={gnoise[19],gnoise[16:15]};
          8 :data_gn[23:21]={gnoise[19],gnoise[16:15]};
          9 :data_gn[20:18]={gnoise[19],gnoise[16:15]};
          10:data_gn[17:15]={gnoise[19],gnoise[16:15]};
          11:data_gn[14:12]={gnoise[19],gnoise[16:15]};
          12:data_gn[11: 9]={gnoise[19],gnoise[16:15]};
          13:data_gn[ 8: 6]={gnoise[19],gnoise[16:15]};
          14:data_gn[ 5: 3]={gnoise[19],gnoise[16:15]};
          15:data_gn[ 2: 0]={gnoise[19],gnoise[16:15]};
        endcase
         case(l)
           0 :  if( data_in1_w [2] ^ data_in_gn1_w [3] )RBER=RBER+1;
           1 :  if( data_in2_w [2] ^ data_in_gn2_w [3] )RBER=RBER+1;
           2 :  if( data_in3_w [2] ^ data_in_gn3_w [3] )RBER=RBER+1;
           3 :  if( data_in4_w [2] ^ data_in_gn4_w [3] )RBER=RBER+1;
           4 :  if( data_in5_w [2] ^ data_in_gn5_w [3] )RBER=RBER+1;
           5 :  if( data_in6_w [2] ^ data_in_gn6_w [3] )RBER=RBER+1;
           6 :  if( data_in7_w [2] ^ data_in_gn7_w [3] )RBER=RBER+1;
           7 :  if( data_in8_w [2] ^ data_in_gn8_w [3] )RBER=RBER+1;
           8 :  if( data_in9_w [2] ^ data_in_gn9_w [3] )RBER=RBER+1;
           9 :  if( data_in10_w[2] ^ data_in_gn10_w[3] )RBER=RBER+1;
           10:  if( data_in11_w[2] ^ data_in_gn11_w[3] )RBER=RBER+1;
           11:  if( data_in12_w[2] ^ data_in_gn12_w[3] )RBER=RBER+1;
           12:  if( data_in13_w[2] ^ data_in_gn13_w[3] )RBER=RBER+1;
           13:  if( data_in14_w[2] ^ data_in_gn14_w[3] )RBER=RBER+1;
           14:  if( data_in15_w[2] ^ data_in_gn15_w[3] )RBER=RBER+1;
           15:  if( data_in16_w[2] ^ data_in_gn16_w[3] )RBER=RBER+1;
        endcase
       end
      start_Ldata=1;
      start_encode=0;
       @(posedge CLK);
       i=0;
    end
    else if(j>block_size+idle_time && j<=2*block_size+idle_time)
    begin
      e=1;
      start_Ldata=1;
      start_encode=1;
      // data_b=packetmem1[i];
         data_b=$random%32768;
        for(l=0;l<16;l=l+1)
       begin
        gnacc=$random%1024;
        gaussnoise=packetmemgn[gnacc];
        gnoise=gaussnoise*var;
        case(l)
          0 :data_gn[47:45]={gnoise[19],gnoise[16:15]};
          1 :data_gn[44:42]={gnoise[19],gnoise[16:15]};
          2 :data_gn[41:39]={gnoise[19],gnoise[16:15]};
          3 :data_gn[38:36]={gnoise[19],gnoise[16:15]};
          4 :data_gn[35:33]={gnoise[19],gnoise[16:15]};
          5 :data_gn[32:30]={gnoise[19],gnoise[16:15]};
          6 :data_gn[29:27]={gnoise[19],gnoise[16:15]};
          7 :data_gn[26:24]={gnoise[19],gnoise[16:15]};
          8 :data_gn[23:21]={gnoise[19],gnoise[16:15]};
          9 :data_gn[20:18]={gnoise[19],gnoise[16:15]};
          10:data_gn[17:15]={gnoise[19],gnoise[16:15]};
          11:data_gn[14:12]={gnoise[19],gnoise[16:15]};
          12:data_gn[11: 9]={gnoise[19],gnoise[16:15]};
          13:data_gn[ 8: 6]={gnoise[19],gnoise[16:15]};
          14:data_gn[ 5: 3]={gnoise[19],gnoise[16:15]};
          15:data_gn[ 2: 0]={gnoise[19],gnoise[16:15]};
        endcase
         case(l)
           0 :  if( data_in1_w [2] ^ data_in_gn1_w [3] )RBER=RBER+1;
           1 :  if( data_in2_w [2] ^ data_in_gn2_w [3] )RBER=RBER+1;
           2 :  if( data_in3_w [2] ^ data_in_gn3_w [3] )RBER=RBER+1;
           3 :  if( data_in4_w [2] ^ data_in_gn4_w [3] )RBER=RBER+1;
           4 :  if( data_in5_w [2] ^ data_in_gn5_w [3] )RBER=RBER+1;
           5 :  if( data_in6_w [2] ^ data_in_gn6_w [3] )RBER=RBER+1;
           6 :  if( data_in7_w [2] ^ data_in_gn7_w [3] )RBER=RBER+1;
           7 :  if( data_in8_w [2] ^ data_in_gn8_w [3] )RBER=RBER+1;
           8 :  if( data_in9_w [2] ^ data_in_gn9_w [3] )RBER=RBER+1;
           9 :  if( data_in10_w[2] ^ data_in_gn10_w[3] )RBER=RBER+1;
           10:  if( data_in11_w[2] ^ data_in_gn11_w[3] )RBER=RBER+1;
           11:  if( data_in12_w[2] ^ data_in_gn12_w[3] )RBER=RBER+1;
           12:  if( data_in13_w[2] ^ data_in_gn13_w[3] )RBER=RBER+1;
           13:  if( data_in14_w[2] ^ data_in_gn14_w[3] )RBER=RBER+1;
           14:  if( data_in15_w[2] ^ data_in_gn15_w[3] )RBER=RBER+1;
           15:  if( data_in16_w[2] ^ data_in_gn16_w[3] )RBER=RBER+1;
        endcase
       end
      packetmem2[i]=data_b;
       @(posedge CLK);
      i=i+1;
    end
    else if(j>2*block_size+idle_time && j<=2*block_size+2*idle_time)
    begin
        for(l=0;l<16;l=l+1)
       begin
        gnacc=$random%1024;
        gaussnoise=packetmemgn[gnacc];
        gnoise=gaussnoise*var;
        case(l)
          0 :data_gn[47:45]={gnoise[19],gnoise[16:15]};
          1 :data_gn[44:42]={gnoise[19],gnoise[16:15]};
          2 :data_gn[41:39]={gnoise[19],gnoise[16:15]};
          3 :data_gn[38:36]={gnoise[19],gnoise[16:15]};
          4 :data_gn[35:33]={gnoise[19],gnoise[16:15]};
          5 :data_gn[32:30]={gnoise[19],gnoise[16:15]};
          6 :data_gn[29:27]={gnoise[19],gnoise[16:15]};
          7 :data_gn[26:24]={gnoise[19],gnoise[16:15]};
          8 :data_gn[23:21]={gnoise[19],gnoise[16:15]};
          9 :data_gn[20:18]={gnoise[19],gnoise[16:15]};
          10:data_gn[17:15]={gnoise[19],gnoise[16:15]};
          11:data_gn[14:12]={gnoise[19],gnoise[16:15]};
          12:data_gn[11: 9]={gnoise[19],gnoise[16:15]};
          13:data_gn[ 8: 6]={gnoise[19],gnoise[16:15]};
          14:data_gn[ 5: 3]={gnoise[19],gnoise[16:15]};
          15:data_gn[ 2: 0]={gnoise[19],gnoise[16:15]};
        endcase
         case(l)
           0 :  if( data_in1_w [2] ^ data_in_gn1_w [3] )RBER=RBER+1;
           1 :  if( data_in2_w [2] ^ data_in_gn2_w [3] )RBER=RBER+1;
           2 :  if( data_in3_w [2] ^ data_in_gn3_w [3] )RBER=RBER+1;
           3 :  if( data_in4_w [2] ^ data_in_gn4_w [3] )RBER=RBER+1;
           4 :  if( data_in5_w [2] ^ data_in_gn5_w [3] )RBER=RBER+1;
           5 :  if( data_in6_w [2] ^ data_in_gn6_w [3] )RBER=RBER+1;
           6 :  if( data_in7_w [2] ^ data_in_gn7_w [3] )RBER=RBER+1;
           7 :  if( data_in8_w [2] ^ data_in_gn8_w [3] )RBER=RBER+1;
           8 :  if( data_in9_w [2] ^ data_in_gn9_w [3] )RBER=RBER+1;
           9 :  if( data_in10_w[2] ^ data_in_gn10_w[3] )RBER=RBER+1;
           10:  if( data_in11_w[2] ^ data_in_gn11_w[3] )RBER=RBER+1;
           11:  if( data_in12_w[2] ^ data_in_gn12_w[3] )RBER=RBER+1;
           12:  if( data_in13_w[2] ^ data_in_gn13_w[3] )RBER=RBER+1;
           13:  if( data_in14_w[2] ^ data_in_gn14_w[3] )RBER=RBER+1;
           14:  if( data_in15_w[2] ^ data_in_gn15_w[3] )RBER=RBER+1;
           15:  if( data_in16_w[2] ^ data_in_gn16_w[3] )RBER=RBER+1;
        endcase
       end
      start_Ldata=1;
      start_encode=0;
       @(posedge CLK);

       i=0;
    end
    else if(j==block_size+idle_time+(m+1)*iteration_time)
    begin
           if(e<2)e=e+1;
           else e=0;
           m=m+1;
        @(posedge CLK);
    end
    else if(j>block_size + idle_time + m*iteration_time && j<=2*block_size + idle_time + m*iteration_time)
    begin
      // data_b=packetmem1[i];
         data_b=$random%32768;
        for(l=0;l<16;l=l+1)
       begin
        gnacc=$random%1024;
        gaussnoise=packetmemgn[gnacc];
        gnoise=gaussnoise*var;
        case(l)
          0 :data_gn[47:45]={gnoise[19],gnoise[16:15]};
          1 :data_gn[44:42]={gnoise[19],gnoise[16:15]};
          2 :data_gn[41:39]={gnoise[19],gnoise[16:15]};
          3 :data_gn[38:36]={gnoise[19],gnoise[16:15]};
          4 :data_gn[35:33]={gnoise[19],gnoise[16:15]};
          5 :data_gn[32:30]={gnoise[19],gnoise[16:15]};
          6 :data_gn[29:27]={gnoise[19],gnoise[16:15]};
          7 :data_gn[26:24]={gnoise[19],gnoise[16:15]};
          8 :data_gn[23:21]={gnoise[19],gnoise[16:15]};
          9 :data_gn[20:18]={gnoise[19],gnoise[16:15]};
          10:data_gn[17:15]={gnoise[19],gnoise[16:15]};
          11:data_gn[14:12]={gnoise[19],gnoise[16:15]};
          12:data_gn[11: 9]={gnoise[19],gnoise[16:15]};
          13:data_gn[ 8: 6]={gnoise[19],gnoise[16:15]};
          14:data_gn[ 5: 3]={gnoise[19],gnoise[16:15]};
          15:data_gn[ 2: 0]={gnoise[19],gnoise[16:15]};
        endcase
         case(l)
           0 :  if( data_in1_w [2] ^ data_in_gn1_w [3] )RBER=RBER+1;
           1 :  if( data_in2_w [2] ^ data_in_gn2_w [3] )RBER=RBER+1;
           2 :  if( data_in3_w [2] ^ data_in_gn3_w [3] )RBER=RBER+1;
           3 :  if( data_in4_w [2] ^ data_in_gn4_w [3] )RBER=RBER+1;
           4 :  if( data_in5_w [2] ^ data_in_gn5_w [3] )RBER=RBER+1;
           5 :  if( data_in6_w [2] ^ data_in_gn6_w [3] )RBER=RBER+1;
           6 :  if( data_in7_w [2] ^ data_in_gn7_w [3] )RBER=RBER+1;
           7 :  if( data_in8_w [2] ^ data_in_gn8_w [3] )RBER=RBER+1;
           8 :  if( data_in9_w [2] ^ data_in_gn9_w [3] )RBER=RBER+1;
           9 :  if( data_in10_w[2] ^ data_in_gn10_w[3] )RBER=RBER+1;
           10:  if( data_in11_w[2] ^ data_in_gn11_w[3] )RBER=RBER+1;
           11:  if( data_in12_w[2] ^ data_in_gn12_w[3] )RBER=RBER+1;
           12:  if( data_in13_w[2] ^ data_in_gn13_w[3] )RBER=RBER+1;
           13:  if( data_in14_w[2] ^ data_in_gn14_w[3] )RBER=RBER+1;
           14:  if( data_in15_w[2] ^ data_in_gn15_w[3] )RBER=RBER+1;
           15:  if( data_in16_w[2] ^ data_in_gn16_w[3] )RBER=RBER+1;
        endcase
       end
      start_Ldata=1;
      start_encode=1;
      if(e==2)packetmem3[i]=data_b;
      else if(e==0)packetmem1[i]=data_b;
      else if(e==1)packetmem2[i]=data_b;
      @(posedge CLK);
      i=i+1;
    end
    else if(j>2*block_size+idle_time+ m*iteration_time && j<=2*block_size+2*idle_time+ m*iteration_time)
    begin
         for(l=0;l<16;l=l+1)
       begin
        gnacc=$random%1024;
        gaussnoise=packetmemgn[gnacc];
        gnoise=gaussnoise*var;
        case(l)
          0 :data_gn[47:45]={gnoise[19],gnoise[16:15]};
          1 :data_gn[44:42]={gnoise[19],gnoise[16:15]};
          2 :data_gn[41:39]={gnoise[19],gnoise[16:15]};
          3 :data_gn[38:36]={gnoise[19],gnoise[16:15]};
          4 :data_gn[35:33]={gnoise[19],gnoise[16:15]};
          5 :data_gn[32:30]={gnoise[19],gnoise[16:15]};
          6 :data_gn[29:27]={gnoise[19],gnoise[16:15]};
          7 :data_gn[26:24]={gnoise[19],gnoise[16:15]};
          8 :data_gn[23:21]={gnoise[19],gnoise[16:15]};
          9 :data_gn[20:18]={gnoise[19],gnoise[16:15]};
          10:data_gn[17:15]={gnoise[19],gnoise[16:15]};
          11:data_gn[14:12]={gnoise[19],gnoise[16:15]};
          12:data_gn[11: 9]={gnoise[19],gnoise[16:15]};
          13:data_gn[ 8: 6]={gnoise[19],gnoise[16:15]};
          14:data_gn[ 5: 3]={gnoise[19],gnoise[16:15]};
          15:data_gn[ 2: 0]={gnoise[19],gnoise[16:15]};
        endcase
          case(l)
           0 :  if( data_in1_w [2] ^ data_in_gn1_w [3] )RBER=RBER+1;
           1 :  if( data_in2_w [2] ^ data_in_gn2_w [3] )RBER=RBER+1;
           2 :  if( data_in3_w [2] ^ data_in_gn3_w [3] )RBER=RBER+1;
           3 :  if( data_in4_w [2] ^ data_in_gn4_w [3] )RBER=RBER+1;
           4 :  if( data_in5_w [2] ^ data_in_gn5_w [3] )RBER=RBER+1;
           5 :  if( data_in6_w [2] ^ data_in_gn6_w [3] )RBER=RBER+1;
           6 :  if( data_in7_w [2] ^ data_in_gn7_w [3] )RBER=RBER+1;
           7 :  if( data_in8_w [2] ^ data_in_gn8_w [3] )RBER=RBER+1;
           8 :  if( data_in9_w [2] ^ data_in_gn9_w [3] )RBER=RBER+1;
           9 :  if( data_in10_w[2] ^ data_in_gn10_w[3] )RBER=RBER+1;
           10:  if( data_in11_w[2] ^ data_in_gn11_w[3] )RBER=RBER+1;
           11:  if( data_in12_w[2] ^ data_in_gn12_w[3] )RBER=RBER+1;
           12:  if( data_in13_w[2] ^ data_in_gn13_w[3] )RBER=RBER+1;
           13:  if( data_in14_w[2] ^ data_in_gn14_w[3] )RBER=RBER+1;
           14:  if( data_in15_w[2] ^ data_in_gn15_w[3] )RBER=RBER+1;
           15:  if( data_in16_w[2] ^ data_in_gn16_w[3] )RBER=RBER+1;
        endcase
       end
      start_Ldata=1;
      start_encode=0;
       @(posedge CLK);
      i=0;
    end
    else
    begin
      start_Ldata=0;
       start_encode=0;
      @(posedge CLK);
      i=0;
    end

         //     Recieve data
    if(hb_vlid)
    begin
         if(e==2)data=packetmem1[k];
         else if(e==0)data=packetmem2[k];
         else if(e==1)data=packetmem3[k];
         data_out=hardbit;
       for (n=0;n<16;n=n+1)
       begin
          if(p<9216)
         if(  hardbit[15-n] != data[n])berr=berr+1;
       end
       k=k+1;
       p=p+16;
    end
    else
    begin
      k=0;
      p=0;
    end
  end
     A=m*9216;
     B=berr;
     C=RBER;
     BER=B/A;
     rbber=C/A;
   $display($time,"%d bit error RBER=%g BER=%g ",berr,rbber,BER);
   $display("%d block Input Fin \n",m);
     RBER=0;
     berr=0;
     m=0;
  end
  $finish;
end
endmodule
