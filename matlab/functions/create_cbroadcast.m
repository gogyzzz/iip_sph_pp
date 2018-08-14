function create_cbroadcast( a0,a1,a2,a3,a4,a5,a6,a7 )

iii = read_cmat(a0);
iid = read_cmat(a1);
idi = read_cmat(a2);
idd= read_cmat(a3);
dii = read_cmat(a4);
did= read_cmat(a5);
ddi= read_cmat(a6);
ddd = read_cmat(a7);

t = '../test/Broadcast.bin';
% 4 4
c = bsxfun(@plus, iid, idd );
f = post_appen(t,'0_cAdd Elements');
write_cmat(f,c);
c = bsxfun(@plus, iid,did );
f = post_appen(t,'1_cAdd Elements');
write_cmat(f,c);
c = bsxfun(@plus, idd,iid );
f = post_appen(t,'2_cAdd Elements');
write_cmat(f,c);
% 4 1qqwwqe
c = bsxfun(@plus, did,iii );
f = post_appen(t,'3_cAdd Elements');
write_cmat(f,c);
c = bsxfun(@times, iid, ddi);
f = post_appen(t,'4_cMul Elements');
write_cmat(f,c);
c = bsxfun(@times, ddd,iii );
f = post_appen(t,'5_cMul Elements');
write_cmat(f,c);
% 1 4
c = bsxfun(@times, idi,did );
f = post_appen(t,'6_cMul Elements');
write_cmat(f,c);
c = bsxfun(@times, dii,idd );
f = post_appen(t,'7_cMul Elements');
write_cmat(f,c);
c = bsxfun(@times, idi, ddd);
f = post_appen(t,'8_cMul Elements');
write_cmat(f,c);
% 1 1
c = bsxfun(@rdivide,dii ,ddi );
f = post_appen(t,'9_cDiv Elements');
write_cmat(f,c);
c = bsxfun(@rdivide, ddi,idi );
f = post_appen(t,'10_cDiv Elements');
write_cmat(f,c);
c = bsxfun(@rdivide, ddi,dii );
f = post_appen(t,'11_cDiv Elements');
write_cmat(f,c);
c = bsxfun(@rdivide, ddi,ddi );
f = post_appen(t,'12_cDiv Elements');
write_cmat(f,c);



end

