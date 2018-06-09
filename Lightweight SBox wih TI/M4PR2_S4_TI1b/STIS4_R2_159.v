module  STIS4_R2_159(in,out);
input[7:0] in;
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
assign term_0=in[0];
assign term_1=in[1];
assign term_2=in[2];
assign term_3=in[0]&in[1];
assign term_4=in[2]&in[3];
assign term_5=in[0]&in[2];
assign term_6=in[1]&in[3];
assign term_7=in[2]&in[4];
assign term_8=in[3]&in[5];
assign term_9=in[1]&in[4];
assign term_10=in[3]&in[6];
assign term_11=in[0]&in[5];
assign term_12=in[2]&in[7];
assign term_13=in[0]&in[6];
assign term_14=in[1]&in[7];
assign out=term_0^term_1^term_2^term_3^term_4^term_5^term_6^term_7^term_8^term_9^term_10^term_11^term_12^term_13^term_14;

endmodule
