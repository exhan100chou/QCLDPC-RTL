# QCLDPC RTL Verilog Channal Codec
for the paper "A novel data packing technique for QC-LDPC decoder architecture applied to NAND flash controller"\
https://ieeexplore.ieee.org/abstract/document/9015393 \

Project: Data Packing QC-LDPC Decoder 
1. (A) Nature order data packing into memory (B) The proposed order data packing 

![alt text](https://github.com/exhan100chou/QCLDPC-RTL/blob/main/photo/Fig1.jpg)

2. QC circularly shifts the data packing in the proposed data packing.
   The red arrow represents that after the cyclic shift of the data pack in line with the check node process. \
   This illustration demonstrates the proposed data packing can avoid data conflict of random QC property    

![alt text](https://github.com/exhan100chou/QCLDPC-RTL/blob/main/photo/Fig2.jpg)

3. The proposed partial parallel architecture 6bit quantization
   
![alt text](https://github.com/exhan100chou/QCLDPC-RTL/blob/main/photo/Fig3.jpg)

4. Finite state machine timesheet 
![alt text](https://github.com/exhan100chou/QCLDPC-RTL/blob/main/photo/Fig4.jpg)

5. \
Testbench 1 of LDPC decoder 48bit input/ 40*8bit output:\
          Input code data from file codeword.txt\
           Test algorithm behavior for each single block\
Testbench 2 of LDPC decoder 48bit input/ 40*8bit output:\
          Gaussian noise gen can run the desired block number\
          Test the BER performance\
