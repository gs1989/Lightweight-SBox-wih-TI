module  STIP8_R3_64876239(in,out);
input[3:0] in;
output[3:0] out;
wire[3:0] out;
assign out[7]=in[6]^in[5];
assign out[6]=in[6]^in[4];
assign out[5]=in[6]^in[3];
assign out[4]=in[6]^in[2];
assign out[3]=in[6]^in[1];
assign out[2]=in[6]^in[0];
assign out[1]=in[6];
assign out[0]=in[7]^in[6];
endmodule
