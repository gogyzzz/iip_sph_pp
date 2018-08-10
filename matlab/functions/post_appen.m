function y = post_appen(filename,appen )
[t,filename] = strtok(filename,'/');
filename = strtok(filename,'/');
c = 'data/';
t = strtok(filename,'.');
c = strcat(c,t);
c =strcat(c,appen);
c =strcat(c,'.bin');
y = c;
end

