/////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    ISPD 2013 sample ///////////////////////////////////
////
//// The netlist from ISPD 2012 contest was placed using Cadence Encounter 
//// Digital Implementation System to produce realistic spef files.
//// All cell sizes were reset to the largest sizes afterwards.
////
////
////////////////////////////////    ISPD 2012 changes ///////////////////////////////////
//// The original netlist was downloaded from http://www.iwls.org/iwls2005/benchmarks.htm
//// and remapped to the ISPD 2012 contest library. Also a few simplification changes 
//// were done on the netlist in order to have :
//// -  no ground/power nets
//// -  no dangling nets,
//// -	no unconnected pin,
//// -	no constants/no assigns
//// -	single clock
//// -	no busses
//// -	no special characters/no escapes/no brackets in names
//// -  fanout optimization 
////
//// The original copyright is included next:
////
////
////  This design was downloaded from http://www.opencores.org
////
////  The design was synthesized with Cadence RTL Compiler in a
////  quick synthesis run.
////
/////////////////////////////////////////////////////////////////////
////                                                             ////
////  USB 1.1 PHY                                                ////
////                                                             ////
////  Author: Rudolf Usselmann                                   ////
////          rudi@asics.ws                                      ////
////                                                             ////
////  Downloaded from: http://www.opencores.org/cores/usb_phy/   ////
////                                                             ////
/////////////////////////////////////////////////////////////////////
////                                                             ////
//// Copyright (C) 2000-2002 Rudolf Usselmann                    ////
////                         www.asics.ws                        ////
////                         rudi@asics.ws                       ////
////                                                             ////
//// This source file may be used and distributed without        ////
//// restriction provided that this copyright statement is not   ////
//// removed from the file and that any derivative work contains ////
//// the original copyright notice and the associated disclaimer.////
////                                                             ////
////     THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY     ////
//// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED   ////
//// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS   ////
//// FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL THE AUTHOR      ////
//// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,         ////
//// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES    ////
//// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE   ////
//// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR        ////
//// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ////
//// LIABILITY, WHETHER IN  CONTRACT, STRICT LIABILITY, OR TORT  ////
//// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT  ////
//// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE         ////
//// POSSIBILITY OF SUCH DAMAGE.                                 ////
////                                                             ////
/////////////////////////////////////////////////////////////////////

// Generated by Cadence RTL Compiler (RC) v05.10-b006_1


module sample (
DataIn_i_0_,
ispd_clk,
DataIn_o_0_,
DataIn_o_1_,
LineState_A_o_,
LineState_B_o_,
rst_cnt_reg_1_o_,
);

// Start PIs
input DataIn_i_0_;
input ispd_clk;

// Start POs
output DataIn_o_0_;
output DataIn_o_1_;

// Start wires
wire DataIn_i_0_;
wire DataIn_o_0_;
wire DataIn_o_1_;
wire LineState_A_o_;
wire LineState_B_o_;
wire ispd_clk;
wire rst_cnt_reg_1_o_;

// Start cells
na02m80 FE_A_0 ( .a(DataIn_i_0_), .b(DataIn_o_0_), .o(LineState_A_o_) );
na02m80 FE_B_0 ( .a(DataIn_o_0_), .b(rst_cnt_reg_1_o_), .o(LineState_B_o_) );
in01f80 FE_C_0 ( .a(LineState_A_o_), .o(DataIn_o_1_) );
ms00f80 rst_cnt_reg_1__u0 ( .ck(ispd_clk), .d(LineState_A_o_), .o(rst_cnt_reg_1_o_) );
ms00f80 rst_cnt_reg_2__u0 ( .ck(ispd_clk), .d(LineState_B_o_), .o(DataIn_o_0_) );

endmodule