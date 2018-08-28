%{
Create output files name <filename>_<function>.bin
ex)
filename = d_4_3_2.bin
output : d_4_3_2_Invert.bin
%}
function y = create_coutput(filename)

A = read_cmat(filename);
total = 0;
[d0 d1 d2] = size(A);

%% Need to be square
if d0 == d1
    
tic;   
    for i = 1:d2;
    outmat(:,:,i) = inv(A(:,:,i));  
    end
    total = total + toc; outname = post_appen(filename,'cInvert');
    write_cmat(outname,outmat);
    
tic;   
    for i = 1:d2;
    outmat(:,:,i) = transpose(A(:,:,i));    
    end
    total = total + toc; outname = post_appen(filename,'cTranspose');
    write_cmat(outname,outmat);

    clearvars outmat
tic;   
    for i = 1:d2;
    t = diag(A(:,:,i));
    outmat(:,i) = t;
    end
    total = total + toc; outname = post_appen(filename,'cDiagonal');
    write_cmat(outname,outmat);

     clearvars outmat
tic;   
    for i = 1:d2;
    outmat(i) = trace(A(:,:,i));
    end
    total = total + toc; outname = post_appen(filename,'cTrace');
    write_cmat(outname,outmat);
    
end


tic;   
outmat = power(A,2);
total = total + toc; outname = post_appen(filename,'Pow cMat');
write_cmat(outname,outmat);

tic;   
outmat = power(A,2 + 2i);
total = total + toc; outname = post_appen(filename,'cPow cMat');
write_cmat(outname,outmat);


tic;  
outmat = sqrt(A);
total = total + toc; outname = post_appen(filename,'cSqrt');
write_cmat(outname,outmat);

tic;   
outmat = round(A);
total = total + toc; outname = post_appen(filename,'cRound');
write_cmat(outname,outmat);

tic;   
outmat = floor(A);
total = total + toc; outname = post_appen(filename,'cFloor');
write_cmat(outname,outmat);

tic;   
outmat = ceil(A);
total = total + toc; outname = post_appen(filename,'cCeil');
write_cmat(outname,outmat);

tic;   
outmat = log(A);
total = total + toc; outname = post_appen(filename,'cLog');
write_cmat(outname,outmat);

tic;   
outmat = log2(A);
total = total + toc; outname = post_appen(filename,'cLog2');
write_cmat(outname,outmat);

tic;   
outmat = log10(A);
total = total + toc; outname = post_appen(filename,'cLog10');
write_cmat(outname,outmat);


tic;   
outmat = exp(A);
total = total + toc; outname = post_appen(filename,'cExp');
write_cmat(outname,outmat);

tic;   
outmat = abs(A);
total = total + toc; outname = post_appen(filename,'cAbs');
write_cmat(outname,outmat);

tic;   
outmat = repmat(A,3,2);
total = total + toc; outname = post_appen(filename,'cRepmat');
write_cmat(outname,outmat);

tic;   
outmat = permute(A,[1 3 2]);
total = total + toc; outname = post_appen(filename,'cPermute132'); % d0 d2 d1
write_cmat(outname,outmat);

tic;   
outmat = permute(A,[2 1 3]);
total = total + toc; outname = post_appen(filename,'cPermute213'); % d1 d0 d2
write_cmat(outname,outmat);

tic;   
outmat = permute(A,[2 3 1]);
total = total + toc; outname = post_appen(filename,'cPermute231'); % d1 d2 d0
write_cmat(outname,outmat);

tic;   
outmat = permute(A,[3, 1 2]);
total = total + toc; outname = post_appen(filename,'cPermute312'); % d2 d0 d1
write_cmat(outname,outmat);

tic;   
outmat = permute(A,[3 2 1]);
total = total + toc; outname = post_appen(filename,'cPermute321'); % d2 d1 d0
write_cmat(outname,outmat);

tic;   
outmat = reshape(A,d0*d1*d2,1,1);
total = total + toc; outname = post_appen(filename,'cReshape'); % d0
write_cmat(outname,outmat);

tic;   
outmat  = A + (10 + 10*1i);
total = total + toc; outname = post_appen(filename,'uAdd');
write_cmat(outname,outmat);

tic;   
outmat  = A + 10 ;
total = total + toc; outname = post_appen(filename,'cAdd');
write_cmat(outname,outmat);

tic;   
outmat = A * 2;
total = total + toc; outname = post_appen(filename,'cScale');
write_cmat(outname,outmat);

tic;   
outmat = A * (2 + 2*1i);
total = total + toc; outname = post_appen(filename,'uScale');
write_cmat(outname,outmat);


clearvars outmat;
tic;   
outmat = sum(A,2);
total = total + toc; outname = post_appen(filename,'Sum cMat'); % output -> d2,1,1
write_cmat(outname,outmat);


clearvars outmat;
tic;   
outmat = sum(A,1);
total = total + toc; outname = post_appen(filename,'Sum cMat 1'); % output -> d2,1,1
write_cmat(outname,outmat);
y = total;
end
