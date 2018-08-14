%{
Create binary FILE named 'filename'
with MATRIX 'mat'
ex)
A = [1 2 ; 3 4]
file data seqeunce 
: 1 3 2 4
%}
function write_mat(filename,mat)
fileID = fopen(filename,'w');
fwrite(fileID, mat, 'double');
fclose(fileID);
disp(filename);
end

