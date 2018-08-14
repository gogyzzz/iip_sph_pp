function create_coutput(filename)

A = read_cmat(filename);

[d0 d1 d2] = size(A);

%% Need to be square
if d0 == d1   
    for i = 1:d2;
    outmat(:,:,i) = inv(A(:,:,i));  
    end
    outname = post_appen(filename,'cInvert');
    write_cmat(outname,outmat);
    
    for i = 1:d2;
    outmat(:,:,i) = transpose(A(:,:,i));    
    end
    outname = post_appen(filename,'cTranspose');
    write_cmat(outname,outmat);

    clearvars outmat
    for i = 1:d2;
    t = diag(A(:,:,i));
    outmat(:,i) = t;
    end
    outname = post_appen(filename,'cDiagonal');
    write_cmat(outname,outmat);

     clearvars outmat
    for i = 1:d2;
    outmat(i) = trace(A(:,:,i));
    end
    outname = post_appen(filename,'cTrace');
    write_cmat(outname,outmat);
    
end


outmat = power(A,2);
outname = post_appen(filename,'cPow');
write_cmat(outname,outmat);

outmat = power(A,2 + 2i);
outname = post_appen(filename,'uPow');
write_cmat(outname,outmat);



outmat = sqrt(A);
outname = post_appen(filename,'cSqrt');
write_cmat(outname,outmat);

outmat = round(A);
outname = post_appen(filename,'cRound');
write_cmat(outname,outmat);

outmat = floor(A);
outname = post_appen(filename,'cFloor');
write_cmat(outname,outmat);

outmat = ceil(A);
outname = post_appen(filename,'cCeil');
write_cmat(outname,outmat);

outmat = log(A);
outname = post_appen(filename,'cLog');
write_cmat(outname,outmat);

outmat = log2(A);
outname = post_appen(filename,'cLog2');
write_cmat(outname,outmat);

outmat = log10(A);
outname = post_appen(filename,'cLog10');
write_cmat(outname,outmat);


outmat = exp(A);
outname = post_appen(filename,'cExp');
write_cmat(outname,outmat);

outmat = abs(A);
outname = post_appen(filename,'cAbs');
write_cmat(outname,outmat);

outmat = repmat(A,3,2);
outname = post_appen(filename,'cRepmat');
write_cmat(outname,outmat);

outmat = repmat(A,3,2);
outname = post_appen(filename,'Repmat'); % d0*3, d1*2 , d2
write_cmat(outname,outmat);

outmat = permute(A,[1 3 2]);
outname = post_appen(filename,'Permute132'); % d0 d2 d1
write_cmat(outname,outmat);

outmat = permute(A,[2 1 3]);
outname = post_appen(filename,'Permute213'); % d1 d0 d2
write_cmat(outname,outmat);

outmat = permute(A,[2 3 1]);
outname = post_appen(filename,'Permute231'); % d1 d2 d0
write_cmat(outname,outmat);

outmat = permute(A,[3, 1 2]);
outname = post_appen(filename,'Permute312'); % d2 d0 d1
write_cmat(outname,outmat);

outmat = permute(A,[3 2 1]);
outname = post_appen(filename,'Permute321'); % d2 d1 d0
write_cmat(outname,outmat);

outmat = reshape(A,d0*d1*d2,1,1);
outname = post_appen(filename,'cReshape'); % d0
write_cmat(outname,outmat);

outmat  = A + (10 + 10*1i);
outname = post_appen(filename,'uAdd');
write_cmat(outname,outmat);

outmat  = A + 10 ;
outname = post_appen(filename,'cAdd');
write_cmat(outname,outmat);

outmat = A * 2;
outname = post_appen(filename,'cScale');
write_cmat(outname,outmat);

outmat = A * (2 + 2*1i);
outname = post_appen(filename,'uScale');
write_cmat(outname,outmat);


for i = 1:d2
outmat(i) = sum(sum(A(:,:,i)));
end
outname = post_appen(filename,'cSum');
write_cmat(outname,outmat);

for i = 1:d2
outmat(i) = sum(sum(A(:,:,i)));
end
outname = post_appen(filename,'cAsum');
write_cmat(outname,outmat);

end
