function write_cmat(filename,mat,d0,d1,d2)
mat_real = real(mat);
mat_imag = imag(mat);
out_mat = zeros(d0*2,d1,d2);

for i = 1:d2
    newrow = 1;
    for row = 1:d0
        out_mat(newrow,:,i) = mat_real(row,:,i);
        out_mat(newrow + 1,:,i) = mat_imag(row,:,i);
        newrow = newrow + 2;
    end
end
fileID = fopen(filename,'w');
fwrite(fileID, d0, 'int');
fwrite(fileID, d1, 'int');
fwrite(fileID, d2, 'int');
fwrite(fileID, out_mat, 'double');
fclose(fileID);
end
