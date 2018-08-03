function [ output_args ] = read_cmat_wip( input_args )
%READ_CMAT_WIP Summary of this function goes here
%   Detailed explanation goes here

 z = zeros(mat_size);
    
    f = fopen(filename);
    %disp(['read from ',filename]);
    %disp(size(z));
    same_real = fread(f, size(z(:)), 'double', 8);
    fseek(f, 8, 'bof');
    same_imag = fread(f, size(z(:)), 'double', 8);
    fclose(f);
    disp(['read ',filename, ' with size ', num2str(size(z(:)))]);
     
    same_z = complex(same_real, same_imag);
    same_z = reshape(same_z, mat_size);
    disp(['convert size ', num2str(size(z(:))), ' to ', num2str(mat_size)]);
%disp(same_z);

end

