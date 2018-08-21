function y=create_broadcast( a0,a1,a2,a3,a4,a5,a6,a7 )

iii = read_mat(a0);
iid = read_mat(a1);
idi = read_mat(a2);
idd= read_mat(a3);
dii = read_mat(a4);
did= read_mat(a5);
ddi= read_mat(a6);
ddd = read_mat(a7);

t = '../test/Broadcast.bin';

total =0;

% 4 4
tic; c = bsxfun(@plus, iid, idd );
total = total + toc; f = post_appen(t,'0_Add Elements');
write_mat(f,c);
tic; c = bsxfun(@plus, iid,did );
total = total + toc; f = post_appen(t,'1_Add Elements');
write_mat(f,c);
tic; c = bsxfun(@plus, idd,iid );
total = total + toc; f = post_appen(t,'2_Add Elements');
write_mat(f,c);
% 4 1
tic; c = bsxfun(@plus, did,iii );
total = total + toc; f = post_appen(t,'3_Add Elements');
write_mat(f,c);
tic; c = bsxfun(@times, iid, ddi);
total = total + toc; f = post_appen(t,'4_Mul Elements');
write_mat(f,c);
tic; c = bsxfun(@times, ddd,iii );
total = total + toc; f = post_appen(t,'5_Mul Elements');
write_mat(f,c);
% 1 4
tic; c = bsxfun(@times, idi,did );
total = total + toc; f = post_appen(t,'6_Mul Elements');
write_mat(f,c);
tic; c = bsxfun(@times, dii,idd );
total = total + toc; f = post_appen(t,'7_Mul Elements');
write_mat(f,c);
tic; c = bsxfun(@times, idi, ddd);
total = total + toc; f = post_appen(t,'8_Mul Elements');
write_mat(f,c);
% 1 1
tic; c = bsxfun(@rdivide,dii ,ddi );
total = total + toc; f = post_appen(t,'9_Div Elements');
write_mat(f,c);
tic; c = bsxfun(@rdivide, ddi,idi );
total = total + toc; f = post_appen(t,'10_Div Elements');
write_mat(f,c);
tic; c = bsxfun(@rdivide, ddi,dii );
total = total + toc; f = post_appen(t,'11_Div Elements');
write_mat(f,c);
tic; c = bsxfun(@rdivide, ddi,ddi );
total = total + toc; f = post_appen(t,'12_Div Elements');
write_mat(f,c);

y = total;

end

