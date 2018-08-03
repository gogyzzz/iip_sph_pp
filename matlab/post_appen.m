function y = post_appen(filename,appen )
c = strtok(filename,'.');
c =strcat(c,appen);
c =strcat(c,'.bin');
y = c;
end

