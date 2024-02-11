//Main module

module MSpVM #(
 parameter INW = 24,
 parameter OUTW = 64,
 parameter M=17,
 parameter N=15,
 localparam LOGN = $clog2(N),
 localparam LOGMN = $clog2(M*N),
 parameter READ=0
 )(
 input clk, reset,
 input [INW-1:0] INPUT_TDATA,
 input INPUT_TVALID,
 input INPUT_TLAST,
 input [LOGN:0] INPUT_TUSER,
 output INPUT_TREADY,
 output [OUTW-1:0] OUTPUT_TDATA,
 output OUTPUT_TVALID,
 input OUTPUT_TREADY
 );


logic  input_loaded, done;
logic [LOGN-1:0] D;
logic [LOGN-1:0] countv;     //-> this is vector_read_addr
logic [INW-1:0] vector_val, vector_val_d;
logic [LOGN-1:0] vector_row;
logic [LOGMN-1:0] matrix_read_addr;
logic [INW-1:0] matrix_data;
logic [OUTW-1:0] data_in_fifo;
logic clear_acc, valid_input;
logic wr_en;
logic [($clog2(M+1))-1:0] capacity;


//
logic incmrow, clrmrow, inccountv, clrcountv;
logic [($clog2(M))-1:0] mrow;
logic [($clog2(M+1))-1:0] mrow1;
logic [($clog2(N+1))-1:0] countv1;

//
logic state, next_state;
//
logic flag, temp1, temp2, temp3, tempdone;

counter1 #(M) mrowcounter(clk, incmrow, clrmrow, mrow);
counter2 #(N) vcounter(clk, inccountv, clrcountv, countv);
counter3 #(M) mrowcounter1(clk, incmrow, clrmrow, mrow1);
counter4 #(N) vcounter1(clk, inccountv, clrcountv, countv1);
 




input_mems #(INW,M,N) input_mem_dut(.clk(clk), .reset(reset), .AXIS_TDATA(INPUT_TDATA), .AXIS_TVALID(INPUT_TVALID),
                                     .AXIS_TLAST(INPUT_TLAST), .AXIS_TUSER(INPUT_TUSER), .AXIS_TREADY(INPUT_TREADY),
                                     .input_loaded(input_loaded), .done(done), .D(D), .vector_read_addr(countv),
                                     .vector_val(vector_val), .vector_row(vector_row), .matrix_read_addr(matrix_read_addr),
                                     .matrix_data(matrix_data) ); 
mac_pipe #(INW,OUTW) mac_pipe_dut (.in0(matrix_data), .in1(vector_val_d), .clk(clk), .reset(reset), .clear_acc(clear_acc), .valid_input(valid_input),
                                    .out(data_in_fifo));

output_fifo #(OUTW,M) fifo_dut (.clk(clk), .reset(reset), .wr_en(wr_en), .capacity(capacity), .data_in(data_in_fifo),
                                     .AXIS_TDATA(OUTPUT_TDATA), .AXIS_TVALID(OUTPUT_TVALID),.AXIS_TREADY(OUTPUT_TREADY));



//
assign matrix_read_addr = mrow*N + vector_row;



always_ff @(posedge clk) begin
    vector_val_d<=vector_val;
    if (reset==1)
        state<=READ;
    else 
        state<=next_state;
end


always_comb begin
        if (input_loaded==0)
            next_state=READ;
        else if (input_loaded==1 && flag==0) 
            next_state=READ;
        else if (flag==1 && valid_input==1) 
            next_state=READ;
        else if (flag==1 && valid_input==0) begin
            if (temp3==0)
                next_state=READ;
            else if (temp3==1 && capacity==0)
                next_state=READ;
            else
                next_state=READ;
        end     
        //else if (tempdone==1)
        else
            next_state=READ;
end


always_ff @(posedge clk) begin
    if (reset)
        temp2<=0;
    else 
        temp2<=temp1;
end

always_ff @(posedge clk) begin
    if (reset)
        temp3<=0;
    else if ((state==READ && flag==1 && valid_input==0) && (temp3==1 && capacity==0))
        temp3<=1;
    else
        temp3<=temp2;
end


always_ff @(posedge clk) begin
    if (reset || tempdone==1 || (flag==1 && valid_input==0 && temp3==1 && capacity>0))
        flag<=0;
    else if ((state==READ && input_loaded==1 && flag==0) && (countv1==D))
        flag<=1;
    else 
        flag<=flag;
end


always_ff @(posedge clk) begin
    if (reset)
            valid_input<=0;
    else if (input_loaded==1 && flag==0 && countv1>0)
        valid_input<=1;
    else
        valid_input<=0;
end

always_ff @(posedge clk) begin
    if (reset==1 || tempdone==1)
        tempdone<=0;
    else if (flag==1 && valid_input==0 && temp3==1 && capacity>0)
        tempdone<=1;
end


assign clear_acc = ((reset==1) || ((flag==1 && valid_input==0) && (temp3==1 && capacity>0)));
assign clrmrow = (reset==1) || (tempdone==1 && mrow1==M);
assign incmrow = ((flag==1 && valid_input==0) && (temp3==1 && capacity>0));
assign clrcountv = (reset==1) || (wr_en==1);
assign inccountv = ((input_loaded==1 && flag==0) && (countv1<D));
assign wr_en = ((flag==1 && valid_input==0) && (temp3==1 && capacity>0));
assign temp1 = (flag==1 && valid_input==1);
assign done = (tempdone==1 && mrow1==M);



endmodule


module counter1(clk, inc, clr, count);
    input clk, clr, inc;
    parameter M = 8;
    localparam LOGSIZE = $clog2(M);
    output logic [LOGSIZE-1:0] count;
    always_ff @(posedge clk) begin
        if (clr == 1)
            count<=0;
        else if (inc)
            count<=count+1;
    end
endmodule

module counter2(clk, inc, clr, count);
    input clk, clr, inc;
    parameter P = 8;
    localparam LOGSIZE = $clog2(P);
    output logic [LOGSIZE-1:0] count;
    always_ff @(posedge clk) begin
        if (clr == 1)
            count<=0;
        else if (inc)
            count<=count+1;
    end
endmodule


module counter3(clk, inc, clr, count);
    input clk, clr, inc;
    parameter Q = 8;
    localparam LOGSIZE = $clog2(Q+1);
    output logic [LOGSIZE-1:0] count;
    always_ff @(posedge clk) begin
        if (clr == 1)
            count<=0;
        else if (inc)
            count<=count+1;
    end
endmodule


module counter4(clk, inc, clr, count);
    input clk, clr, inc;
    parameter R = 8;
    localparam LOGSIZE = $clog2(R+1);
    output logic [LOGSIZE-1:0] count;
    always_ff @(posedge clk) begin
        if (clr == 1)
            count<=0;
        else if (inc)
            count<=count+1;
    end
endmodule





//--------------------




//PART 3 INPUT MEMORY
module memory #(parameter WIDTH=16, SIZE=32, localparam LOGSIZE=$clog2(SIZE)) 
(input [WIDTH-1:0] data_in, output logic [WIDTH-1:0] data_out, input [LOGSIZE-1:0] addr, input clk, wr_en);

	logic [SIZE-1:0] [WIDTH-1:0] mem;
	always_ff @(posedge clk) begin
		data_out <= mem[addr];
		if (wr_en)
			mem[addr] <= data_in;
	end
endmodule


module matrix_ctr(clk, inc_m, clear_m, m_ctr);
	input clk, clear_m, inc_m;
	parameter M = 3, N= 3;
	localparam LOGSIZE = $clog2(M*N);
	output logic [LOGSIZE-1:0] m_ctr;
	always_ff @(posedge clk) begin
		if (clear_m == 1)
			m_ctr<=0;
		else if (inc_m)
			m_ctr<=m_ctr+1;
	end
endmodule


module vector_ctr(clk, inc_v, clear_v, v_ctr);
	input clk, clear_v, inc_v;
	parameter X =3;
	localparam LOGSIZE = $clog2(X);
	output logic [LOGSIZE-1:0] v_ctr;
	always_ff @(posedge clk) begin
		if (clear_v == 1)
			v_ctr<=0;
		else if (inc_v)
			v_ctr<=v_ctr+1;
	end
endmodule



module input_mems #(parameter INW=8, M=5, N=5, localparam LOGN = $clog2(N), LOGMN = $clog2(M*N))
(input clk, reset,
input [INW-1:0] AXIS_TDATA,
input AXIS_TVALID,
input AXIS_TLAST,
input [LOGN:0] AXIS_TUSER,
output logic AXIS_TREADY,
output logic input_loaded,
input done,
output logic [LOGN-1:0] D,
input [LOGN-1:0] vector_read_addr,
output logic [INW-1:0] vector_val,
output logic [LOGN-1:0] vector_row,
input [LOGMN-1:0] matrix_read_addr,
output logic [INW-1:0] matrix_data); 


logic new_matrix;
assign new_matrix = AXIS_TUSER[0];

enum {initital_state, store_state, loaded_state}state, next_state;

always_ff @(posedge clk) begin
	if (reset==1)
		state<=initital_state;
	else 
		state<=next_state;
end

localparam LOGMN1 = $clog2(M*N+1);
localparam LOGN1 = $clog2(N+1);

logic wr_en_m, wr_en_v;
logic [LOGMN1-1:0] maddr;
logic [LOGN1-1:0] vaddr;
logic [LOGMN-1:0] mataddr;
logic [LOGN-1:0] vecaddr;
logic [LOGMN-1:0] maddr1;
logic [LOGN-1:0] vaddr1; 
logic [(INW+LOGN)-1:0] vin, vout;
assign vin = {AXIS_TUSER[LOGN:1], AXIS_TDATA};

always_comb begin
	if (input_loaded==1) begin
		mataddr = matrix_read_addr;
		vecaddr = vector_read_addr;
	end
	else begin
		mataddr = maddr1;
		vecaddr = vaddr1;
	end
end

memory #(INW, M*N) mat(AXIS_TDATA, matrix_data, mataddr, clk, wr_en_m);
memory #(INW+LOGN, N) vec(vin, vout, vecaddr, clk, wr_en_v);
assign vector_row = vout [(INW+LOGN)-1:INW];
assign vector_val = vout [INW-1:0];
	

logic incrmaddr, clearmaddr;
logic incrvaddr, clearvaddr;
matrix_ctr #(M, N) matrix_ctrInst(clk, incrmaddr, clearmaddr, maddr1);
vector_ctr #(N) vector_ctrInst(clk, incrvaddr, clearvaddr, vaddr1);
vector_ctr #(M*N+1) vector_ctrInst2(clk, incrmaddr, clearmaddr, maddr);
vector_ctr #(N+1) vector_ctrInst1(clk, incrvaddr, clearvaddr, vaddr);


always_ff @(posedge clk) begin
	if (reset==1)
		D<=0;
	else if (state==store_state && AXIS_TLAST==1)
		D<=vaddr+1;
end


always_comb begin
case(state)
	initital_state:begin
		next_state = store_state;
	end	
	store_state: begin
		if (AXIS_TLAST == 1)
			next_state = loaded_state;
		else 
			next_state = store_state;
	end
	loaded_state: begin
		if (done == 1)
			next_state = store_state;
		else
			next_state = loaded_state;
	end
	default:begin
		next_state = initital_state;
	end	
endcase		
end
	
assign clearmaddr = (state == initital_state) || (state == loaded_state && done==1);
assign clearvaddr = (state == initital_state) || (state == loaded_state && done==1);
assign AXIS_TREADY = (state == store_state) ;
assign wr_en_m = (state==store_state && AXIS_TLAST==0 && AXIS_TVALID==1 && ((maddr==0 && new_matrix==1 && vaddr==0) || (maddr>0 && maddr<M*N)));
assign incrmaddr = (state==store_state && AXIS_TLAST==0 && AXIS_TVALID==1 && ((maddr==0 && new_matrix==1 && vaddr==0) || (maddr>0 && maddr<M*N)));
assign wr_en_v = (state == store_state && AXIS_TLAST==0 && AXIS_TVALID==1 && ((maddr==0 && vaddr==0 && new_matrix==0) || (vaddr>0) || (maddr>M*N-1) )) || (state == store_state && AXIS_TLAST==1);
assign incrvaddr = (state == store_state && AXIS_TLAST==0 && AXIS_TVALID==1 && ((maddr==0 && vaddr==0 && new_matrix==0) || (vaddr>0) || (maddr>M*N-1) )) || (state == store_state && AXIS_TLAST==1);
assign input_loaded = (state == loaded_state && done==0);

endmodule



//---------------


//PART1 MAC Unit

module mac_pipe #(
 parameter INW = 16,
 parameter OUTW = 64
 )(
 input signed [INW-1:0] in0, in1,
 output logic signed [OUTW-1:0] out,
 input clk, reset, clear_acc, valid_input
 );

 logic signed [2*INW-1:0] p;
 logic signed [2*INW-1:0] pr;
 logic signed [OUTW-1:0] f;
 logic en;
 
 always_comb begin
	p= in0*in1;
	f = pr+ out;
 end

 
 always_ff @(posedge clk) begin
	if (reset == 1 || clear_acc==1)	
		pr <= 0;
	else
		pr <= p;
 end 

 always_ff @(posedge clk) begin
 	if (valid_input==1)
 		en<=1;
 	else
 		en<=0;
 end


 always_ff @(posedge clk) begin
	if (reset == 1 || clear_acc==1)	
		out <= 0;
	else if (en == 1)
		out <= f;
 end 
 
 
 
endmodule
	


//--------------


//PART2 FIFO Unit

module output_fifo (clk, reset, data_in, wr_en, capacity, AXIS_TDATA, AXIS_TVALID, AXIS_TREADY);
	parameter OUTW=28, DEPTH=8;
	localparam LOGDEPTH=$clog2(DEPTH);
	
	input clk, reset, wr_en, AXIS_TREADY;
	input [OUTW-1:0] data_in;
	output logic [OUTW-1:0] AXIS_TDATA;
	output logic AXIS_TVALID;
	output logic [($clog2(DEPTH+1))-1:0] capacity;
	
	logic clearhead, incrhead, cleartail, clearcap, incrcap, decrcap;
	logic incrtail;
	logic [LOGDEPTH-1:0] counthead;
	logic [LOGDEPTH-1:0] counttail, addr;
	logic rd_en;
	
	memory_dual_port #(OUTW, DEPTH) memdual(data_in, AXIS_TDATA, counthead, addr, clk, wr_en);
	
	head #(DEPTH) counterhead (clk, clearhead, incrhead, counthead);
	tail #(DEPTH) countertail (clk, cleartail, incrtail, counttail);
	capacity #(DEPTH) countercapacity (clk, clearcap, incrcap, decrcap, capacity);
	
	
	parameter START = 0, OPERATION = 1;
	logic state, nextState;
	
	always_ff @(posedge clk) begin
		if (reset == 1)
			state <= START;
		else
			state <= nextState;
	end
	
	always_comb begin
		if (state == START)
			nextState = OPERATION;
		else if (state == OPERATION)
			nextState = OPERATION;
		else 
			nextState = START;
	end
	
	always_comb begin
		if (rd_en)
			if (counttail==DEPTH-1)
				addr = 0;
			else 
				addr = counttail+1;
		else 
			addr = counttail;
	end
	
	
	assign AXIS_TVALID = ((state==OPERATION) && (capacity<DEPTH));
	assign rd_en = (AXIS_TREADY==1 && AXIS_TVALID==1);
	assign clearhead = (state==START);
	assign cleartail = (state==START);
	assign clearcap = (state==START);
	assign incrhead = (state==OPERATION && wr_en==1);
	assign incrtail = (state==OPERATION && rd_en==1);
	assign incrcap = (state==OPERATION && rd_en==1 && wr_en==0);
	assign decrcap = (state==OPERATION && wr_en==1 && rd_en==0);
	
	
endmodule 


module memory_dual_port #(
 parameter WIDTH=16, SIZE=64,
 localparam LOGSIZE=$clog2(SIZE)
 )(
 input [WIDTH-1:0] data_in,
 output logic [WIDTH-1:0] data_out,
 input [LOGSIZE-1:0] write_addr, read_addr,
 input clk, wr_en
 );

 logic [SIZE-1:0][WIDTH-1:0] mem;

 always_ff @(posedge clk) begin
 data_out <= mem[read_addr];
 if (wr_en) begin
 mem[write_addr] <= data_in;
 if (read_addr == write_addr)
 data_out <= data_in;
 end
 end
endmodule



module head(clk, clearhead, incrhead, counthead);
	parameter H = 8;
	localparam LOGSIZE = $clog2(H);
	
	input clk, clearhead, incrhead;
	output logic [LOGSIZE-1:0] counthead;
	
	always_ff @(posedge clk) begin
		if (clearhead == 1)
			counthead<=0;
		else if (incrhead==1) begin
			if (counthead==H-1)
				counthead<=0;
			else
				counthead<=counthead+1;
		end
	end
	
endmodule
	

module tail(clk, cleartail, incrtail, counttail);
	parameter T =8;
	localparam LOGSIZE = $clog2(T);
	
	input clk, cleartail, incrtail;
	output logic [LOGSIZE-1:0] counttail;
	
	always_ff @(posedge clk) begin
		if (cleartail == 1)
			counttail<=0;
		else if (incrtail==1) begin
			if (counttail==T-1)
				counttail<=0;
			else
				counttail<=counttail+1;
		end
	end

endmodule	
	

module capacity(clk, clearcap, incrcap, decrcap, countcap);
	parameter C =8;
	localparam LOGSIZE = $clog2(C);
	
	input clk, clearcap, incrcap, decrcap;
	output logic [($clog2(C+1))-1:0] countcap;
	
	always_ff @(posedge clk) begin
		if (clearcap == 1)
			countcap<=C;
		else if (incrcap==1)
			countcap<=countcap+1;
		else if (decrcap==1)
			countcap<=countcap-1;
	end
	
endmodule	


