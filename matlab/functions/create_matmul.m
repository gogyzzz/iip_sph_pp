function y = create_matmul(f_a,f_b,f_c)
A = read_mat(f_a);
B = read_mat(f_b);

[d0 d1 d2] = size(A);
tic;
for i = 1:d2
C(:,:,i) = A(:,:,i) * B(:,:,i);
end
y = toc;
write_mat(f_c,C);

end

