module  STIS8_R3_34433786(in,out);
input[15:0] in;
output out;
wire out;
wire term_0;
wire term_1;
wire term_2;
wire term_3;
wire term_4;
wire term_5;
wire term_6;
wire term_7;
wire term_8;
wire term_9;
wire term_10;
wire term_11;
wire term_12;
wire term_13;
wire term_14;
wire term_15;
wire term_16;
wire term_17;
wire term_18;
wire term_19;
wire term_20;
wire term_21;
wire term_22;
wire term_23;
wire term_24;
wire term_25;
wire term_26;
wire term_27;
wire term_28;
wire term_29;
wire term_30;
wire term_31;
wire term_32;
wire term_33;
wire term_34;
wire term_35;
wire term_36;
wire term_37;
wire term_38;
assign term_0=in[1];
assign term_1=in[3];
assign term_2=in[6];
assign term_3=in[0]&in[1];
assign term_4=in[1]&in[2];
assign term_5=in[2]&in[3];
assign term_6=in[5]&in[6];
assign term_7=in[0]&in[2];
assign term_8=in[1]&in[3];
assign term_9=in[6]&in[8];
assign term_10=in[7]&in[9];
assign term_11=in[1]&in[4];
assign term_12=in[7]&in[10];
assign term_13=in[0]&in[4];
assign term_14=in[1]&in[5];
assign term_15=in[4]&in[8];
assign term_16=in[5]&in[9];
assign term_17=in[2]&in[7];
assign term_18=in[4]&in[9];
assign term_19=in[0]&in[6];
assign term_20=in[1]&in[7];
assign term_21=in[2]&in[8];
assign term_22=in[3]&in[9];
assign term_23=in[1]&in[8];
assign term_24=in[2]&in[9];
assign term_25=in[3]&in[10];
assign term_26=in[6]&in[13];
assign term_27=in[0]&in[9];
assign term_28=in[1]&in[10];
assign term_29=in[2]&in[11];
assign term_30=in[5]&in[14];
assign term_31=in[0]&in[10];
assign term_32=in[1]&in[11];
assign term_33=in[1]&in[12];
assign term_34=in[0]&in[12];
assign term_35=in[1]&in[13];
assign term_36=in[2]&in[15];
assign term_37=in[0]&in[14];
assign term_38=in[1]&in[15];
assign out=term_0^term_1^term_2^term_3^term_4^term_5^term_6^term_7^term_8^term_9^term_10^term_11^term_12^term_13^term_14^term_15^term_16^term_17^term_18^term_19^term_20^term_21^term_22^term_23^term_24^term_25^term_26^term_27^term_28^term_29^term_30^term_31^term_32^term_33^term_34^term_35^term_36^term_37^term_38;

endmodule
