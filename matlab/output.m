
create_output('../test_data/d_1024_1024_4.bin');
%create_output('base/d_1_1_1.bin');
%create_matmul('base/d_4_1024_1.bin','base/d_1024_6_1.bin','data/d_4_6_1_cMatmul.bin');
create_matmul('../test_data/d_1_1024_4.bin','../test_data/d_1024_1024_4.bin','../test_ans/d_1_1024_4_Gemm.bin');
create_broadcast('../test_data/d_1_1_1.bin','../test_data/d_1_1_4.bin','../test_data/d_1_1024_1.bin','../test_data/d_1_1024_4.bin','../test_data/d_1024_1_1.bin','../test_data/d_1024_1_4.bin','../test_data/d_1024_1024_1.bin','../test_data/d_1024_1024_4.bin');

create_coutput('../test_data/c_1024_1024_4.bin');
%create_coutput('base/c_1_1_1.bin');
%create_cmatmul('base/c_4_1024_1.bin','base/c_1024_6_1.bin','data/c_4_6_1_cMatmul.bin');
create_cmatmul('../test_data/c_1_1024_4.bin','../test_data/c_1024_1024_4.bin','../test_ans/c_1_1024_4_cGemm.bin');
create_cbroadcast('../test_data/c_1_1_1.bin','../test_data/c_1_1_4.bin','../test_data/c_1_1024_1.bin','../test_data/c_1_1024_4.bin','../test_data/c_1024_1_1.bin','../test_data/c_1024_1_4.bin','../test_data/c_1024_1024_1.bin','../test_data/c_1024_1024_4.bin');

%{
create_output('../test_data/d_4_4_2.bin');
%create_output('base/d_1_1_1.bin');
%create_matmul('base/d_2_4_1.bin','base/d_4_6_1.bin','data/d_2_6_1_cMatmul.bin');
create_matmul('../test_data/d_1_4_2.bin','../test_data/d_4_4_2.bin','../test_ans/d_1_4_2_Gemm.bin');
create_broadcast('../test_data/d_1_1_1.bin','../test_data/d_1_1_2.bin','../test_data/d_1_4_1.bin','../test_data/d_1_4_2.bin','../test_data/d_4_1_1.bin','../test_data/d_4_1_2.bin','../test_data/d_4_4_1.bin','../test_data/d_4_4_2.bin');

create_coutput('../test_data/c_4_4_2.bin');
%create_coutput('base/c_1_1_1.bin');
%create_cmatmul('base/c_2_4_1.bin','base/c_4_6_1.bin','data/c_2_6_1_cMatmul.bin');
create_cmatmul('../test_data/c_1_4_2.bin','../test_data/c_4_4_2.bin','../test_ans/c_1_4_2_cGemm.bin');
create_cbroadcast('../test_data/c_1_1_1.bin','../test_data/c_1_1_2.bin','../test_data/c_1_4_1.bin','../test_data/c_1_4_2.bin','../test_data/c_4_1_1.bin','../test_data/c_4_1_2.bin','../test_data/c_4_4_1.bin','../test_data/c_4_4_2.bin');
%}