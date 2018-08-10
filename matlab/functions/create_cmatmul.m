function create_cmatmul(f_a,f_b,f_c)
A = read_cmat(f_a);
B = read_cmat(f_b);

[d0 d1 d2] = size(A);

for i = 1:d2
C(:,:,i) = A(:,:,i) * B(:,:,i);
end
write_cmat(f_c,C);

end

