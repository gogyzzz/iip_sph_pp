function y = create_mat(d0,d1,d2)
col = d0;
row = d1;
batch = d2;
A = rand(col,row,batch); 
A = A*100;
y=A;
f = 'base/d_';
t = num2str(col);
f=strcat(f,t);
f=strcat(f,'_');
t = num2str(row);
f=strcat(f,t);
f=strcat(f,'_');
t = num2str(batch);
f=strcat(f,t);
f=strcat(f,'.bin');
disp(f);
write_mat(f,A);

end