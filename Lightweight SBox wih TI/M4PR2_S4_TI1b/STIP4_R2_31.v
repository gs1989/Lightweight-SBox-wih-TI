module  STIP4_R2_31(in,out);
input[3:0] in;
output[3:0] out;
wire[3:0] out;
assign out[3]=in[2]^in[1];
assign out[2]=in[2]^in[0];
assign out[1]=in[2];
assign out[0]=in[3]^in[2];
endmodule
