function y = time_matmul( A,B )
[d0 d1 d2] = size(A);
tic;
for i = 1:d2
C(:,:,i) = A(:,:,i) * B(:,:,i);
end
y = toc; +10 -10
end

