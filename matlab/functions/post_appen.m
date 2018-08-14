% ex) filename = d_4_3_2.bin, appen = hello
% y = data/d_4_3_2_hello.bin
function y = post_appen(filename,appen )
[t1,filename] = strtok(filename,'/');
[t1,filename] = strtok(filename,'/');
filename = strtok(filename,'/');
c = '../test_ans/';
t = strtok(filename,'.');
c = strcat(c,t);
c =strcat(c,'_');
c =strcat(c,appen);
c =strcat(c,'.bin');
y = c;
end

