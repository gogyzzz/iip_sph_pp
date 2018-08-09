function y = post_appen(filename,appen )
c = 'data/';
t = strtok(filename,'.');
c = strcat(c,t);
c =strcat(c,appen);
c =strcat(c,'.bin');
y = c;
end

