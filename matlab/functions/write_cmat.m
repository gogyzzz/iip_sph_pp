%{
Create binary FILE named 'filename'
with complex MATRIX 'mat'
ex)
A = [1 + 2i, 3+ 4i; 5+ 6i, 7+8i;]
file data sequence 
: 1 2 5 6 3 4 7 8
%}
function write_cmat(filename,mat)

mat_real = real(mat);
mat_imag = imag(mat);
[d0,d1,d2] = size(mat);

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
fwrite(fileID, out_mat, 'double');
fclose(fileID);
disp(filename);
end
