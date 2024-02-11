// To use this testbench:

// Compile it and your accompanying design with:
//   vlog -64 +floatparameters +acc MSpVM_tb.sv [add your other .sv files to simulate here]
//   vsim -64 -c MSpVM_tb -sv_seed random 
//      [options]:
//       - If you want to run in GUI mode, remove -c
//       - from the command line you can set parameters in your top-level module by adding:         -g <PARAMNAME>=<VALUE>
//              
//         For example, if you want to set: M=5, N=7, INPUT_TVALID_PROB=0.4, and OUTPUT_TREADY_PROB=0.6, then you could run the following:
//      
// vsim -64 -c MSpVM_tb -g M=5 -g N=7 -g INPUT_TVALID_PROB=0.403 -g OUTPUT_TREADY_PROB=0.6


// Import the C function (in test_helper.c) that generates a random sparse input vector as an array.
import "DPI-C" function void genVector(output int vec_vals[], output int vec_rows[], input int D, input int INW);

// Import the C function (in test_helper.c) that computes the output vector given a matrix and sparse input vector.
import "DPI-C" function void calcOutput(input int matrix[], input int M, input int N, input int vec_vals[], input int vec_rows[], input int D, output longint outvec[], input int OUTW);


// A class to hold one instance of test data and its expected output. When we call 
// .randomize() on an object of this class, it will randomly generate values for the 
// matrix, the number of non-zeros in the vector, and the value of new_matrix.
// This class also has two functions force_new_matrix() and allow_old_matrix().
// These functions will generate random input values and calculate the expected
// output vector. The difference between these is that force_new_matrix() will
// always set new_matrix=1 while allow_old_matrix() will randomly pick whether to
// use the old matrix or create a new one.
class testdata #(parameter INW=8, parameter OUTW=32, parameter M=18, parameter N=15);
    
    // our MxN matrix
    rand int matrix[M*N];   
    
    // This constraint will ensure that the matrix values fit within INW bits.
    constraint m {
        foreach(matrix[i]) {
            matrix[i] < (1<<(INW-1));
            matrix[i] >= -1*(1<<(INW-1));                        
        }
    }

    // vector_vals and vector_rows represent our sparse vector in compressed format.
    // Note that we allocate N entries in this array, but our vector 
    // will use D entries (where D is randomly picked to be 1 <= D <= N).
    // The random vector will be generated using a C function genVector(...) called
    // via DPI. (See input_mems_tb.sv from Part 3 for more explanation)
    int vector_vals[N];      
    int vector_rows[N]; 
    
    // The number of non-zero entries in our sparse vector
    rand logic [$clog2(N)-1:0] D;

    // This constraint means when we randomize an object of this class, it will
    // choose the value of D to be between 1 and N, with equal probability of 
    // each number.
    constraint c {D dist {[1:N] := 1};}
    
    // new_matrix==0 means this test will use the previous matrix; new_matrix==1 means this
    // test will load a new input matrix.
    rand logic new_matrix;    

    // The expected output vector of this MSpVM calculation
    longint output_vector[M];

    // This will randomize the input matrix and sparse vector and compute the 
    // expected result. This function will force new_matrix to be 1.     
    function void force_new_matrix();
        // randomize matrix, D, and new_matrix
        assert(this.randomize());  

        // force new_matrix to 1
        this.new_matrix = 1;       
        
        // generate sparse input vector
        genVector(this.vector_vals, this.vector_rows, this.D, INW);
        
        // calculate the expected result and store in this.output_vector
        calcOutput(this.matrix, M, N, this.vector_vals, this.vector_rows, this.D, this.output_vector, OUTW);
    endfunction

    // This will randomly choose new_matrix. If it is 1, it will generate a new
    // random matrix. If it is 0, it will copy the old matrix from olddata.
    // Then it will generate a new sparse input vector and compute the
    // expected output.
    function void allow_old_matrix(testdata #(INW, OUTW, M, N) olddata);
        // randomize matrix, D, and new_matrix
        assert(this.randomize());

        // if the random new_matrix is 0, then we will copy the previous matrix
        // from old_data
        if (this.new_matrix == 0) 
            this.matrix = olddata.matrix;

        // generate sparse input vector
        genVector(this.vector_vals, this.vector_rows, this.D, INW);

        // calculate the expected result and store in this.output_vector     
        calcOutput(this.matrix, M, N, this.vector_vals, this.vector_rows, this.D, this.output_vector, OUTW);

    endfunction

endclass


module MSpVM_tb();

    parameter TESTS = 10000; // the number of MVMs to test
    parameter INW  = 24;       // the number of bits in the input matrix and vector entries. 
    parameter OUTW = 64;       // the number of bits in the output entries

    // We require 2 <= INW < 32 and 4 <= OUTW <= 64.
    // OUTW must also be large enough to prevent the accumulator from overflowing.
    // (If the accumulator overflows, the testbench will warn you when it computes
    // the expected result.)

    parameter M=17;            // the number of rows in the matrix
    parameter N=15;            // the number of cols in the matrix + rows in vector         

    // The probability that the testbench asserts INPUT_TVALID=1 and OUTPUT_TREADY
    // on any given cycle.
    // You can adjust these values to simulate different scenarios.
    // Valid values for these parameters are 0.001 to 1. 
    // If a value is set to 0, then it will be randomized when you start
    // your simulation.
    parameter real INPUT_TVALID_PROB = 0.2;
    parameter real OUTPUT_TREADY_PROB = 0.8;


    localparam LOGN      = $clog2(N);
    localparam LOGMN     = $clog2(M*N);

    logic clk, reset;
    
    // Signals for the DUT's AXI-Stream input interface
    logic signed [INW-1:0] INPUT_TDATA;
    logic INPUT_TVALID;
    logic INPUT_TLAST;
    logic [LOGN:0] INPUT_TUSER; // INPUT_TUSER[LOGN:1] is the encoding of each vector entry.
                                  // INPUT_TUSER[0] is the new_matrix signal
    logic INPUT_TREADY;   

    // Signals for the DUT's AXI-Stream output interface
    logic signed [OUTW-1:0] OUTPUT_TDATA;
    logic                   OUTPUT_TVALID;
    logic                   OUTPUT_TREADY;

    initial clk=0;
    always #5 clk = ~clk;

    // Instantiate the DUT
    MSpVM #(INW, OUTW, M, N) dut(clk, reset, INPUT_TDATA, INPUT_TVALID, INPUT_TLAST, INPUT_TUSER, INPUT_TREADY, OUTPUT_TDATA, OUTPUT_TVALID, OUTPUT_TREADY);

    // This is an array of "testdata" objects. (See definition of the testdata class above.)
    // It will will hold all of our test data. Each "testdata" object holds data for one MSpVM operation test case.
    // It holds a set of inputs (matrix, vector, new_matrix, etc.) and the corresponding expected output. 
    testdata #(INW, OUTW, M, N) td[TESTS];

    // generate random bits to use when randomizing INPUT_TVALID and OUTPUT_TREADY
    logic rb0, rb1;
    logic [9:0] randomNum;
    logic [9:0] tvalid_prob, tready_prob;

    initial begin
        if (INPUT_TVALID_PROB >= 0.001)
            tvalid_prob = (1024*INPUT_TVALID_PROB-1);
        else
            tvalid_prob = ($urandom % 1024);

        if (OUTPUT_TREADY_PROB >= 0.001)
            tready_prob = (1024*OUTPUT_TREADY_PROB-1);
        else
            tready_prob = ($urandom % 1024);

        $display("--------------------------------------------------------");
        $display("Starting top-level simulation: %d tests", TESTS);
        $display("Matrix dimensions: %d x %d", M, N);
        $display("Input %d bits\nOutput: %d bits", INW, OUTW);

        $display("INPUT_TVALID_PROB = %1.3f", real'(tvalid_prob+1)/1024);
        $display("OUTPUT_TREADY_PROB = %1.3f", real'(tready_prob+1)/1024);            
        $display("--------------------------------------------------------");
    end

    // Every clock cycle, randomly generate rb0 and rb1 for the INPUT_TVALID and
    // OUTPUT_TREADY signals, respectively
    always begin
        @(posedge clk);
        #1;
        randomNum = $urandom;
        rb0 = (randomNum <= tvalid_prob);
        randomNum = $urandom;
        rb1 = (randomNum <= tready_prob);
    end

    // Logic to keep track of where we are in the test data.
    // which_test keeps track of which of the TESTS test cases we are operating on
    // which_element keeps track of which matrix/vector value within this test case we
    // are sending.
    logic [31:0] which_test, which_element; 
    initial which_test=0; 
    initial which_element=0;
    always @(posedge clk) begin
        if (INPUT_TVALID && INPUT_TREADY) begin
            if (which_element == M*N+td[which_test].D-1) begin // if we just finished loading this test...
                which_test <= #1 which_test+1;   // increment to next test

                // if we are not at the last test input:
                if (which_test < TESTS-1) begin
                    // if the next test has a new_matrix, set the counter back to 0
                    if (td[which_test+1].new_matrix == 1) begin
                        which_element <= #1 0;
                    end
                    else begin // if doesn't have a new_matrix, set the counter to the vector location
                        which_element <= #1 M*N;
                    end
                end

            end
            else begin   // or if we are in the middle of a test, just increment
                which_element <= #1 which_element+1;
            end        
        end
    end

    // Logic to set INPUT_TVALID based on random value rb0.
    always @* begin
        // If we haven't finished all of our test inputs and the random bit rb0 is 1,
        if ((which_test < TESTS) && (rb0==1'b1))
            INPUT_TVALID=1;
        else
            INPUT_TVALID=0;
    end

    // Logic to set the value of INPUT_TDATA, INPUT_TUSER, and INPUT_TLAST, based
    // on the random input value and the which_element and which_test variables
    always @* begin
        INPUT_TUSER = 'x;
        INPUT_TDATA = 'x;
        INPUT_TLAST = 'x;
        if (INPUT_TVALID == 1) begin
            if (which_element < M*N) begin  // we are loading the matrix
                INPUT_TDATA = td[which_test].matrix[which_element];
                INPUT_TUSER[0] = td[which_test].new_matrix;   
                INPUT_TLAST = 0;
            end
            else begin // we are loading the vector
                INPUT_TDATA = td[which_test].vector_vals[which_element-M*N];
                INPUT_TUSER = {td[which_test].vector_rows[which_element-M*N], td[which_test].new_matrix};
                if (which_element == M*N + td[which_test].D-1)
                    INPUT_TLAST = 1;
                else
                    INPUT_TLAST = 0;                    
            end            
        end
    end    

    // generate our test input data and expected output data
    initial begin        
        td[0]=new();
        td[0].force_new_matrix(); // the first test needs a new_matrix

        for (int i=1; i<TESTS; i++) begin
            td[i]=new();
            td[i].allow_old_matrix(td[i-1]);      
        end
    end

    // Logic to set OUTPUT_TREADY based on random value rb1
    logic [31:0] which_test_out, which_element_out; 
    always @* begin
        if ((which_test_out < TESTS) && (rb1==1'b1))
            OUTPUT_TREADY = 1;
        else
            OUTPUT_TREADY = 0;
    end

    
    integer errors = 0;
    initial which_test_out = 0;
    initial which_element_out = 0;

    integer cycle_count=0;

    // Logic to check the outputs and keep track of which output test you are checking
    always @(posedge clk) begin
        if (OUTPUT_TVALID && OUTPUT_TREADY) begin 
            if (OUTPUT_TDATA !== td[which_test_out].output_vector[which_element_out]) begin
                $display($time,,"ERROR: Test %d, y[%d] = %d; expected value = %d", which_test_out, which_element_out, OUTPUT_TDATA, td[which_test_out].output_vector[which_element_out]);        
                errors = errors+1;
            end
            if (which_element_out == M-1) begin
                which_element_out = 0;
                which_test_out = which_test_out+1;
            end 
            else begin
                which_element_out = which_element_out+1;
            end
        end
    end

    // Logic to count cycles, used in our throughput testing
    always @(posedge clk) begin

        // reset the cycle_counter on the first element of the first input
        if (INPUT_TVALID && INPUT_TREADY && (which_test==0) && (which_element==0))
            cycle_count <= 0;    
        else
            cycle_count <= cycle_count+1;

    end

    // Logic to assert reset at the beginning, then wait until all tests are done,
    // print the results, and then quit the simulation.
    initial begin
        reset = 1;        
        @(posedge clk); #1; reset = 0; 

        wait(which_test_out == TESTS);
        $display("Simulated %d tests, with a total of %d output values. Detected %d errors.", TESTS, TESTS*N, errors);
        $display("Your system computed %d MSpVMs in %d cycles", TESTS, cycle_count);
        #1;
        $finish;
    end
  

endmodule


