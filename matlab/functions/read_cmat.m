function y = read_cmat(f)
%READ_CMAT_WIP Summary of this function goes here
%   Detailed explanation goes here

if exist(f,'file') 
    bin = fopen(f,'r');
    msg = f;
    msg = strcat(msg,' : reading');
    disp(msg);
      
    [t,f] =strtok(f,'/');
    f = strtok(f,'/');
    c = strtok(f,'.');
    [c,rc] = strtok(c,'_');
    [d0,rc] = strtok(rc,'_');
    [d1,rc] = strtok(rc,'_');
    [d2,rc] = strtok(rc,'_');
    d0 = str2num(d0);
    d1 = str2num(d1);
    d2 = str2num(d2);    

    %mat = zeros(d0,d1,d2);
    temp = fread(bin,'double');
    fclose(bin);
    temp = reshape(temp,d0*2,d1,d2);
   
    mat = zeros(d0,d1,d2);
    for i = 1:d0
        mat(i,:,:) = temp(2*i-1,:,:) + temp(2*i,:,:)*1i;        
    end
    
    
   y =mat;
else
    msg = f;
    msg = strcat(msg,' does not exist');
    disp(msg);
end
end

