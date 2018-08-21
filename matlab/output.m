N = 256;

total = 0.;
d_iii = '../test_data/d_1_1_1.bin';
d_iid = '../test_data/d_1_1_4.bin';
f_d_idi = '../test_data/d_1_%d_1.bin';
f_d_idd = '../test_data/d_1_%d_4.bin';
f_d_dii = '../test_data/d_%d_1_1.bin';
f_d_did = '../test_data/d_%d_1_4.bin';
f_d_ddi = '../test_data/d_%d_%d_1.bin';
f_d_ddd = '../test_data/d_%d_%d_4.bin';

f_d_mul = '../test_ans/d_%d_%d_4_Gemm.bin';

d_idi = sprintf(f_d_idi,N);
d_idd = sprintf(f_d_idd,N);
d_dii = sprintf(f_d_dii,N);
d_did = sprintf(f_d_did,N);
d_ddi = sprintf(f_d_ddi,N,N);
d_ddd = sprintf(f_d_ddd,N,N);

d_mul = sprintf(f_d_mul,N,N);

total = total + create_output(d_ddd);
total = total + create_matmul(d_ddd,d_ddd,d_mul);
 total = total + create_broadcast(d_iii,d_iid,d_idi,d_idd,d_dii,d_did,d_ddi,d_ddd);

c_iii = '../test_data/c_1_1_1.bin';
c_iid = '../test_data/c_1_1_4.bin';
f_c_idi = '../test_data/c_1_%d_1.bin';
f_c_idd = '../test_data/c_1_%d_4.bin';
f_c_dii = '../test_data/c_%d_1_1.bin';
f_c_did = '../test_data/c_%d_1_4.bin';
f_c_ddi = '../test_data/c_%d_%d_1.bin';
f_c_ddd = '../test_data/c_%d_%d_4.bin';

f_c_mul = '../test_ans/c_%d_%d_4_cGemm.bin';

c_idi = sprintf(f_c_idi,N);
c_idd = sprintf(f_c_idd,N);
c_dii = sprintf(f_c_dii,N);
c_did = sprintf(f_c_did,N);
c_ddi = sprintf(f_c_ddi,N,N);
c_ddd = sprintf(f_c_ddd,N,N);

c_mul = sprintf(f_c_mul,N,N);

 total = total + create_coutput(c_ddd);
 total = total + create_cmatmul(c_ddd,c_ddd,c_mul);
 total = total + create_cbroadcast(c_iii,c_iid,c_idi,c_idd,c_dii,c_did,c_ddi,c_ddd);
 disp(total*1000);
