function write_mat(filename,mat)
fileID = fopen(filename,'w');
fwrite(fileID, mat, 'double');
fclose(fileID);
disp(filename);
end

