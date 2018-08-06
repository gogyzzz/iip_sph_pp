
%% READ test_file.txt
test = fopen('test.txt','r');


%% FOR each test file
while ~feof(test)
    %% RUN functions
        test_op  = fscanf(test,'%s',1);
        datatype = fscanf(test,'%c',1);
        d0       = fscanf(test,'%d',1);
        d1       = fscanf(test,'%d',1);
        d2       = fscanf(test,'%d',1);
        %filename = fgetl(test);
        infile = fscanf(test,'%s',1);
        outfile = fscanf(test,'%s',1); 
        A = read_mat(filename); 
        
        
        
        %{
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
          %}
end
fclose(test);

