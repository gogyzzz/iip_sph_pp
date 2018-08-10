function create_output(filename)

A = read_mat(filename);

[d0 d1 d2] = size(A);

%% Need to be square
if d0 == d1
    
    for i = 1:d2;
    outmat(:,:,i) = inv(A(:,:,i));
    end
    outname = post_appen(filename,'_inv');
    write_mat(outname,outmat);
    
    for i = 1:d2;
    outmat(:,:,i) = transpose(A(:,:,i));  
    end
      outname = post_appen(filename,'_trans');
    write_mat(outname,outmat);


    for i = 1:d2;
    outmat(:,:,i) = A(:,:,i)  * A(:,:,i);
    end
    outname = post_appen(filename,'_matmulself');
    write_mat(outname,outmat);

    for i = 1:d2;
    outmat(:,i) = diag(A(:,:,i));
    end
    outname = post_appen(filename,'_diag');
    write_mat(outname,outmat);

    for i = 1:d2;
    outmat(i) = trace(A(:,:,i));
    end
    outname = post_appen(filename,'_trace');
    write_mat(outname,outmat);

end

outmat = power(A,2);
outname = post_appen(filename,'_pow2');
write_mat(outname,outmat);

outmat = sqrt(A);
outname = post_appen(filename,'_sqrt');
write_mat(outname,outmat);

outmat = round(A);
outname = post_appen(filename,'_round');
write_mat(outname,outmat);

outmat = floor(A);
outname = post_appen(filename,'_floor');
write_mat(outname,outmat);

outmat = ceil(A);
outname = post_appen(filename,'_ceil');
write_mat(outname,outmat);

outmat = log(A);
outname = post_appen(filename,'_log');
write_mat(outname,outmat);

outmat = log2(A);
outname = post_appen(filename,'_log2');
write_mat(outname,outmat);

outmat = log10(A);
outname = post_appen(filename,'_log10');
write_mat(outname,outmat);

outmat = exp(A);
outname = post_appen(filename,'_exp');
write_mat(outname,outmat);

outmat = abs(A);
outname = post_appen(filename,'_abs');
write_mat(outname,outmat);

outmat = repmat(A,3,2);
outname = post_appen(filename,'_repmat32');
write_mat(outname,outmat);

outmat = permute(A,[2 3 1]);
outname = post_appen(filename,'_permute231');
write_mat(outname,outmat);

outmat = permute(A,[3, 1 2]);
outname = post_appen(filename,'_permute321');
write_mat(outname,outmat);

outmat = reshape(A,d0*d1*d2,1,1);
outname = post_appen(filename,'_reshaped0');
write_mat(outname,outmat);

outmat = A + 10;
outname = post_appen(filename,'_addp10');
write_mat(outname,outmat);

outmat = A - 10;
outname = post_appen(filename,'_addm10');
write_mat(outname,outmat);

outmat = A * 2;
outname = post_appen(filename,'_scal2');
write_mat(outname,outmat);

outmat = A  * 0.5;
outname = post_appen(filename,'_scalpoint5');
write_mat(outname,outmat);

for i = 1:d2
outmat(i) = sum(sum(A(:,:,i)));
end
outname = post_appen(filename,'_sum');
write_mat(outname,outmat);




end
