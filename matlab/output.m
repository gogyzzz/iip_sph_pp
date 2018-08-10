%{
create_output('base/d_1024_1024_8.bin');
create_output('base/d_1_1_1.bin');
create_matmul('base/d_512_1024_1.bin','base/d_1024_2048_1.bin','data/d_514_2048_1_mul.bin');
create_matmul('base/d_128_256_4.bin','base/d_256_512_4.bin','data/d_128_512_4_mul.bin');
create_broadcast('base/d_1_1_1.bin','base/d_1_1_4.bin','base/d_1_128_1.bin','base/d_1_128_4.bin','base/d_128_1_1.bin','base/d_128_1_4.bin','base/d_128_128_1.bin','base/d_128_128_4.bin');


create_coutput('base/c_1024_1024_8.bin');
create_coutput('base/c_1_1_1.bin');
create_cmatmul('base/c_512_1024_1.bin','base/c_1024_2048_1.bin','data/c_514_2048_1_mul.bin');
create_cmatmul('base/c_128_256_4.bin','base/c_256_512_4.bin','data/c_128_512_4_mul.bin');
create_cbroadcast('base/c_1_1_1.bin','base/c_1_1_4.bin','base/c_1_128_1.bin','base/c_1_128_4.bin','base/c_128_1_1.bin','base/c_128_1_4.bin','base/c_128_128_1.bin','base/c_128_128_4.bin');
%}

create_output('base/d_4_4_2.bin');
create_output('base/d_1_1_1.bin');
create_matmul('base/d_2_4_1.bin','base/d_4_6_1.bin','data/d_2_6_1_mul.bin');
create_matmul('base/d_1_2_2.bin','base/d_2_3_2.bin','data/d_1_3_2_mul.bin');
create_broadcast('base/d_1_1_1.bin','base/d_1_1_2.bin','base/d_1_4_1.bin','base/d_1_4_2.bin','base/d_4_1_1.bin','base/d_4_1_2.bin','base/d_4_4_1.bin','base/d_4_4_2.bin');

create_coutput('base/c_4_4_2.bin');
create_coutput('base/c_1_1_1.bin');
create_cmatmul('base/c_2_4_1.bin','base/c_4_6_1.bin','data/c_2_6_1_mul.bin');
create_cmatmul('base/c_1_2_2.bin','base/c_2_3_2.bin','data/c_1_3_2_mul.bin');
create_cbroadcast('base/c_1_1_1.bin','base/c_1_1_2.bin','base/c_1_4_1.bin','base/c_1_4_2.bin','base/c_4_1_1.bin','base/c_4_1_2.bin','base/c_4_4_1.bin','base/c_4_4_2.bin');

