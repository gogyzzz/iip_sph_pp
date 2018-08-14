function create_broadcast( a0,a1,a2,a3,a4,a5,a6,a7 )

iii = read_mat(a0);
iid = read_mat(a1);
idi = read_mat(a2);
idd= read_mat(a3);
dii = read_mat(a4);
did= read_mat(a5);
ddi= read_mat(a6);
ddd = read_mat(a7);

t = '../test/Broadcast.bin';

% 4 4
c = bsxfun(@plus, iid, idd );
f = post_appen(t,'0_Add Elements');
write_mat(f,c);
c = bsxfun(@plus, iid,did );
f = post_appen(t,'1_Add Elements');
write_mat(f,c);
c = bsxfun(@plus, idd,iid );
f = post_appen(t,'2_Add Elements');
write_mat(f,c);
% 4 1
c = bsxfun(@plus, did,iii );
f = post_appen(t,'3_Add Elements');
write_mat(f,c);
c = bsxfun(@times, iid, ddi);
f = post_appen(t,'4_Mul Elements');
write_mat(f,c);
c = bsxfun(@times, ddd,iii )
f = post_appen(t,'5_Mul Elements');
write_mat(f,c);
% 1 4
c = bsxfun(@times, idi,did )
f = post_appen(t,'6_Mul Elements');
write_mat(f,c);
c = bsxfun(@times, dii,idd )
f = post_appen(t,'7_Mul Elements');
write_mat(f,c);
c = bsxfun(@times, idi, ddd)
f = post_appen(t,'8_Mul Elements');
write_mat(f,c);
% 1 1
c = bsxfun(@rdivide,dii ,ddi );
f = post_appen(t,'9_Div Elements');
write_mat(f,c);
c = bsxfun(@rdivide, ddi,idi );
f = post_appen(t,'10_Div Elements');
write_mat(f,c);
c = bsxfun(@rdivide, ddi,dii );
f = post_appen(t,'11_Div Elements');
write_mat(f,c);
c = bsxfun(@rdivide, ddi,ddi );
f = post_appen(t,'12_Div Elements');
write_mat(f,c);



end

