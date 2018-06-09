module  STIP8_R3_34441753(in,out);
input[3:0] in;
output[3:0] out;
wire[3:0] out;
assign out[7]=in[5];
assign out[6]=in[4];
assign out[5]=in[3];
assign out[4]=in[2];
assign out[3]=in[1];
assign out[2]=in[0];
assign out[1]=in[7];
assign out[0]=in[6];
endmodule
