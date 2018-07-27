* Test circuit version
.subckt ADDER3 A0 A1 A2 B0 B1 B2 CIN S0 S1 S2 COUT VDD GND
x_fa10 A0 B0 GND sum10 cout10 VDD GND FA1
x_fa11 A1 B1 cout10 sum11 cout11 VDD GND FA2 
x_fa12 A2 B2 cout11 sum12 cout12 VDD GND FA3
x_fa20 A0 B0 VDD sum20 cout20 VDD GND FA1
x_fa21 A1 B1 cout20 sum21 cout21 VDD GND FA2
x_fa22 A2 B2 cout21 sum22 cout22 VDD GND FA3
x_mux0 sum10 sum20 CIN S0 VDD GND MUX21
x_mux1 sum11 sum21 CIN S1 VDD GND MUX21
x_mux2 sum12 sum22 CIN S2 VDD GND MUX21
x_muxout cout12 cout22 CIN COUT VDD GND MUX21
.ends and

.subckt FA1 a b c_in sum c_out VDD GND
x_nand1 a b out1 VDD GND NAND_4X
x_nand2 out1 a out2 VDD GND NAND_4X
x_nand3 out1 b out3 VDD GND NAND_4X
x_nand4 out3 out2 out4 VDD GND NAND_4X
x_nand5 out4 c_in out5 VDD GND NAND_4X
x_nand6	out5 out4 out6 VDD GND NAND_3X
x_nand7	out5 c_in out7 VDD GND NAND_3X
x_nand8 out7 out6 sum VDD GND NAND_3X
x_nand9	out5 out1 c_out VDD GND NAND_3X
.ends and

.subckt FA2 a b c_in sum c_out VDD GND
x_nand1 a b out1 VDD GND NAND_3X
x_nand2 out1 a out2 VDD GND NAND_3X
x_nand3 out1 b out3 VDD GND NAND_3X
x_nand4 out3 out2 out4 VDD GND NAND_3X
x_nand5 out4 c_in out5 VDD GND NAND_3X
x_nand6	out5 out4 out6 VDD GND NAND_2X
x_nand7	out5 c_in out7 VDD GND NAND_2X
x_nand8 out7 out6 sum VDD GND NAND_2X
x_nand9	out5 out1 c_out VDD GND NAND_2X
.ends and

.subckt FA3 a b c_in sum c_out VDD GND
x_nand1 a b out1 VDD GND NAND_2X
x_nand2 out1 a out2 VDD GND NAND_2X
x_nand3 out1 b out3 VDD GND NAND_2X
x_nand4 out3 out2 out4 VDD GND NAND_2X
x_nand5 out4 c_in out5 VDD GND NAND_2X
x_nand6	out5 out4 out6 VDD GND NAND_1X
x_nand7	out5 c_in out7 VDD GND NAND_1X
x_nand8 out7 out6 sum VDD GND NAND_1X
x_nand9	out5 out1 c_out VDD GND NAND_1X
.ends and

.subckt MUX21 in1 in2 ctrl out VDD GND
x_mux1_inv	ctrl			INV_out VDD GND INV_2X
x_mux2_nand	in1		INV_out N2_out 	VDD GND NAND_2X
x_mux3_nand	ctrl	in2		N3_out 	VDD GND NAND_2X
x_mux4_nand	N2_out	N3_out	out 	VDD GND NAND_2X
.ends and

.subckt NAND_1X IN1 IN2 OUT VDD GND
mp1 OUT IN1 VDD VDD pch w=0.56u l=0.18u
C0000 OUT 0 0.000816P
C0001 VDD 0 0.000816P
C0002 IN1 0 0.000272P
mp2 OUT IN2 VDD VDD pch w=0.56u l=0.18u
C0003 OUT 0 0.000816P
C0004 VDD 0 0.000816P
C0005 IN2 0 0.000272P
mn1 NET IN1 GND GND nch w=0.56u l=0.18u
C0006 NET 0 0.000816P
C0007 GND 0 0.000816P
C0008 IN1 0 0.000272P
mn2 OUT IN2 NET GND nch w=0.56u l=0.18u
C0009 OUT 0 0.000816P
C0010 NET 0 0.000816P
C0011 IN2 0 0.000272P
.ends and

.subckt NAND_2X IN1 IN2 OUT VDD GND
mp1 OUT IN1 VDD VDD pch w=0.67u l=0.18u
C0012 OUT 0 0.000977P
C0013 VDD 0 0.000977P
C0014 IN1 0 0.000326P
mp2 OUT IN2 VDD VDD pch w=0.67u l=0.18u
C0015 OUT 0 0.000977P
C0016 VDD 0 0.000977P
C0017 IN2 0 0.000326P
mn1 NET IN1 GND GND nch w=0.67u l=0.18u
C0018 NET 0 0.000977P
C0019 GND 0 0.000977P
C0020 IN1 0 0.000326P
mn2 OUT IN2 NET GND nch w=0.67u l=0.18u
C0021 OUT 0 0.000977P
C0022 NET 0 0.000977P
C0023 IN2 0 0.000326P
.ends and

.subckt NAND_3X IN1 IN2 OUT VDD GND
mp1 OUT IN1 VDD VDD pch w=0.80u l=0.18u
C0024 OUT 0 0.001166P
C0025 VDD 0 0.001166P
C0026 IN1 0 0.000389P
mp2 OUT IN2 VDD VDD pch w=0.80u l=0.18u
C0027 OUT 0 0.001166P
C0028 VDD 0 0.001166P
C0029 IN2 0 0.000389P
mn1 NET IN1 GND GND nch w=0.80u l=0.18u
C0030 NET 0 0.001166P
C0031 GND 0 0.001166P
C0032 IN1 0 0.000389P
mn2 OUT IN2 NET GND nch w=0.80u l=0.18u
C0033 OUT 0 0.001166P
C0034 NET 0 0.001166P
C0035 IN2 0 0.000389P
.ends and


.subckt NAND_4X IN1 IN2 OUT VDD GND
mp1 OUT IN1 VDD VDD pch w=1.0u l=0.18u
C0036 OUT 0 0.001458P
C0037 VDD 0 0.001458P
C0038 IN1 0 0.000486P
mp2 OUT IN2 VDD VDD pch w=1.0u l=0.18u
C0039 OUT 0 0.001458P
C0040 VDD 0 0.001458P
C0041 IN2 0 0.000486P
mn1 NET IN1 GND GND nch w=1.0u l=0.18u
C0042 NET 0 0.001458P
C0043 GND 0 0.001458P
C0044 IN1 0 0.000486P
mn2 OUT IN2 NET GND nch w=1.0u l=0.18u
C0045 OUT 0 0.001458P
C0046 NET 0 0.001458P
C0047 IN2 0 0.000486P
.ends and

.subckt INV_1X IN OUT VDD GND
mp1 OUT IN VDD VDD pch w=0.56u l=0.18u
C0048 OUT 0 0.000816P
C0049 VDD 0 0.000816P
C0050 IN 0 0.000272P
mn1 OUT IN GND GND nch w=0.56u l=0.18u
C0051 OUT 0 0.000816P
C0052 GND 0 0.000816P
C0053 IN 0 0.000272P
.ends and


.subckt INV_2X IN OUT VDD GND
mp1 OUT IN VDD VDD pch w=0.67u l=0.18u
C0054 OUT 0 0.000977P
C0055 VDD 0 0.000977P
C0056 IN 0 0.000326P
mn1 OUT IN GND GND nch w=0.67u l=0.18u
C0057 OUT 0 0.000977P
C0058 GND 0 0.000977P
C0059 IN 0 0.000326P
.ends and

.subckt INV_3X IN OUT VDD GND
mp1 OUT IN VDD VDD pch w=0.80u l=0.18u
C0060 OUT 0 0.001166P
C0061 VDD 0 0.001166P
C0062 IN 0 0.000389P
mn1 OUT IN GND GND nch w=0.80u l=0.18u
C0063 OUT 0 0.001166P
C0064 GND 0 0.001166P
C0065 IN 0 0.000389P
.ends and

.subckt INV_4X IN OUT VDD GND
mp1 OUT IN VDD VDD pch w=1.0u l=0.18u
C0066 OUT 0 0.001458P
C0067 VDD 0 0.001458P
C0068 IN 0 0.000486P
mn1 OUT IN GND GND nch w=1.0u l=0.18u
C0069 OUT 0 0.001458P
C0070 GND 0 0.001458P
C0071 IN 0 0.000486P
.ends


VA0 A0 0 PWL ( 5n 0v, 5.05n 1.2v )
VA1 A1 0 PWL ( 5n 0v, 5.05n 1.2v)
VA2 A2 0 PWL ( 5n 0v, 5.05n 1.2v )
VB0 B0 0 PWL ( 5n 0v, 30n 0v )
VB1 B1 0 PWL ( 5n 0v, 30n 0v   )
VB2 B2 0 PWL ( 5n 0v, 30n 0v   )
VCIN CIN 0 PWL ( 5n 0v, 5.05n 1.2v   )
.OPTIONS LIST NODE POST
.TRAN 20P 3100N



.PRINT V(S0) V(S1) V(S2) V(COUT)
.PROB V(A0) V(A1) V(A2) V(B0) V(B1) V(B2) V(CIN) V(S0) V(S1) V(S2) V(COUT) V(A_S0) V(A_S1) V(A_S2) V(A_COUT)
*-----------------------------------------
XADDER A0 A1 A2 B0 B1 B2 CIN S0 S1 S2 COUT VDD GND ADDER3
*-----------------------------------------
VVDD VDD 0 1.2
VGND GND 0 0



CLOADCOUT COUT 0 50f
CLOADS0 S0 0 50f
CLOADS1 S1 0 50f
CLOADS2 S2 0 50f
.MODEL pch PMOS LEVEL=1
.MODEL nch NMOS LEVEL=1
.END