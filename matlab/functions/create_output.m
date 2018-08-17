%{
Create output files name <filename>_<function>.bin
ex)
filename = d_4_3_2.bin
output : d_4_3_2_Invert.bin
%}
function create_output(filename)

A = read_mat(filename);

[d0 d1 d2] = size(A);

%% Need to be square
if d0 == d1
    
    for i = 1:d2;
    outmat(:,:,i) = inv(A(:,:,i));
    end
    outname = post_appen(filename,'Invert');
    write_mat(outname,outmat);
    
    for i = 1:d2;
    outmat(:,:,i) = transpose(A(:,:,i));  
    end
      outname = post_appen(filename,'Transpose'); % d1 d0 d2
    write_mat(outname,outmat);

    clearvars outmat
    for i = 1:d2;
    t = diag(A(:,:,i));
    outmat(:,i) = t;
    end
    outname = post_appen(filename,'Diagonal'); % d0 1 d2
    write_mat(outname,outmat);

    clearvars outmat
    
    for i = 1:d2;
    outmat(i) = trace(A(:,:,i));
    end
    outname = post_appen(filename,'Trace'); % 1 1 d2
    write_mat(outname,outmat);

end

outmat = power(A,2);
outname = post_appen(filename,'Pow Mat');
write_mat(outname,outmat);

outmat = sqrt(abs(A));
outname = post_appen(filename,'Sqrt');
write_mat(outname,outmat);

outmat = round(A);
outname = post_appen(filename,'Round');
write_mat(outname,outmat);

outmat = floor(A);
outname = post_appen(filename,'Floor');
write_mat(outname,outmat);

outmat = ceil(A);
outname = post_appen(filename,'Ceil');
write_mat(outname,outmat);

outmat = log(abs(A));
outname = post_appen(filename,'Log');
write_mat(outname,outmat);

outmat = log2(abs(A));
outname = post_appen(filename,'Log2');
write_mat(outname,outmat);

outmat = log10(abs(A));
outname = post_appen(filename,'Log10');
write_mat(outname,outmat);

outmat = exp(abs(A));
outname = post_appen(filename,'Exp');
write_mat(outname,outmat);

outmat = abs(A);
outname = post_appen(filename,'Abs');
write_mat(outname,outmat);

outmat = repmat(A,3,2);
outname = post_appen(filename,'Repmat'); % d0*3, d1*2 , d2
write_mat(outname,outmat);

outmat = permute(A,[1 3 2]);
outname = post_appen(filename,'Permute132'); % d0 d2 d1
write_mat(outname,outmat);

outmat = permute(A,[2 1 3]);
outname = post_appen(filename,'Permute213'); % d1 d0 d2
write_mat(outname,outmat);

outmat = permute(A,[2 3 1]);
outname = post_appen(filename,'Permute231'); % d1 d2 d0
write_mat(outname,outmat);

outmat = permute(A,[3, 1 2]);
outname = post_appen(filename,'Permute312'); % d2 d0 d1
write_mat(outname,outmat);

outmat = permute(A,[3 2 1]);
outname = post_appen(filename,'Permute321'); % d2 d1 d0
write_mat(outname,outmat);

outmat = reshape(A,d0*d1*d2,1,1);
outname = post_appen(filename,'Reshape'); %d0*d1*d2 1 1
write_mat(outname,outmat);

outmat = A + 10;
outname = post_appen(filename,'Add');
write_mat(outname,outmat);

outmat = A * 2;
outname = post_appen(filename,'Scale');
write_mat(outname,outmat);

clearvars outmat;
outmat = sum(A,1);
outname = post_appen(filename,'Sum Mat 1'); % output -> d0,1,d2
write_mat(outname,outmat);


clearvars outmat;
outmat = sum(A,2);
outname = post_appen(filename,'Sum Mat 0'); % output -> d1,d0,d2
write_mat(outname,outmat);

end
