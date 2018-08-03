function write_mat(filename,mat,d0,d1,d2)
fileID = fopen(filename,'w');
fwrite(fileID, mat, 'double');
fclose(fileID);
end

