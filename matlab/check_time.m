%% REAL TYPE

A0 =read_mat('base/d_1024_1024_8.bin');
t0 = time_output(A0)

A0 =read_mat('base/d_1_1_1.bin');
t1 = time_output(A0)

A0 = read_mat('base/d_512_1024_1.bin');
A1 = read_mat('base/d_1024_2048_1.bin');
t2 = time_matmul(A0,A1)

A0 = read_mat('base/d_128_256_4.bin');
A1 = read_mat('base/d_256_512_4.bin');
t3 = time_matmul(A0,A1)

A0 = read_mat('base/d_1_1_1.bin');
A1 = read_mat('base/d_1_1_4.bin');
A2 = read_mat('base/d_1_128_1.bin');
A3 = read_mat('base/d_1_128_4.bin');
A4 = read_mat('base/d_128_1_1.bin');
A5 = read_mat('base/d_128_1_4.bin');
A6 = read_mat('base/d_128_128_1.bin');
A7 = read_mat('base/d_128_128_4.bin');
t4 = time_broadcast(A0,A1,A2,A3,A4,A5,A6,A7)

total = t0 + t1 + t2 + t3 + t4

%% COMPLEX TYPE
A0 =read_cmat('base/c_1024_1024_8.bin');
t0 = time_output(A0)

A0 =read_cmat('base/c_1_1_1.bin');
t1 = time_output(A0)

A0 = read_cmat('base/c_512_1024_1.bin');
A1 = read_cmat('base/c_1024_2048_1.bin');
t2 = time_matmul(A0,A1)

A0 = read_cmat('base/c_128_256_4.bin');
A1 = read_cmat('base/c_256_512_4.bin');
t3 = time_matmul(A0,A1)

A0 = read_cmat('base/c_1_1_1.bin');
A1 = read_cmat('base/c_1_1_4.bin');
A2 = read_cmat('base/c_1_128_1.bin');
A3 = read_cmat('base/c_1_128_4.bin');
A4 = read_cmat('base/c_128_1_1.bin');
A5 = read_cmat('base/c_128_1_4.bin');
A6 = read_cmat('base/c_128_128_1.bin');
A7 = read_cmat('base/c_128_128_4.bin');
t4 = time_broadcast(A0,A1,A2,A3,A4,A5,A6,A7)

total = t0 + t1 + t2 + t3 + t4