
%% READ test_file.txt
test = fopen('test.txt','r');


%% FOR each test file
while ~feof(test)
    %% RUN functions
        filename = fgetl(test);
        A = read_mat(filename);     
    %% SAVE output files    
        %% IsIt SQAURE ?        
        if size(A,1) == size(A,2)
            %% Invert
            B = zeros(size(A));
            for i = 1:size(A,3);
                B(:,:,i) = inv(A(:,:,i))   ;         
            end
            name = post_appen(filename,'_inv');
            write_mat(name,B,size(B,1),size(B,2),size(B,3));            
            %%          
          
         end
          %% Tranpose
            B = permute(A,[2 1 3]);
            name = post_appen(filename,'_tran');
            write_mat(name,B,size(B,1),size(B,2),size(B,3));
            
          %% Matmul
          % How would I do?
          
          %% Permute
          % Tese all cases? 
          
end
fclose(test);

