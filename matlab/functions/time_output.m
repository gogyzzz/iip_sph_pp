function y = time_output( A )
tic
[d0 d1 d2] = size(A);


%% Need to be square
if d0 == d1
    
    for i = 1:d2;
    outmat(:,:,i) = inv(A(:,:,i));
    end
    
    for i = 1:d2;
    outmat(:,:,i) = transpose(A(:,:,i));  
    end


    for i = 1:d2;
    outmat(:,:,i) = A(:,:,i)  * A(:,:,i);
    end

    for i = 1:d2;
    outmat(:,i) = diag(A(:,:,i));
    end
    for i = 1:d2;
    outmat(i) = trace(A(:,:,i));
    end
end

outmat = power(A,2);
outmat = sqrt(A);

outmat = round(A);

outmat = floor(A);

outmat = ceil(A);
outmat = log(A);

outmat = log2(A);
outmat = log10(A);

outmat = exp(A);
outmat = abs(A);

outmat = repmat(A,3,2);

outmat = permute(A,[2 3 1]);

outmat = permute(A,[3, 1 2]);

outmat = reshape(A,d0*d1*d2,1,1);

outmat = A + 10;

outmat = A - 10;

outmat = A * 2;

outmat = A  * 0.5;

for i = 1:d2
outmat(i) = sum(sum(A(:,:,i)));
end

y = toc;

end

