function y = create_cbroadcast( a0,a1,a2,a3,a4,a5,a6,a7 )
total = 0;
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
tic; c = bsxfun(@plus, iid, idd );
total = total + toc; f = post_appen(t,'0_cAdd Elements');
write_cmat(f,c);
tic; c = bsxfun(@plus, iid,did );
total = total + toc; f = post_appen(t,'1_cAdd Elements');
write_cmat(f,c);
tic; c = bsxfun(@plus, idd,iid );
total = total + toc; f = post_appen(t,'2_cAdd Elements');
write_cmat(f,c);
% 4 1qqwwqe
tic; c = bsxfun(@plus, did,iii );
total = total + toc; f = post_appen(t,'3_cAdd Elements');
write_cmat(f,c);
tic; c = bsxfun(@times, iid, ddi);
total = total + toc; f = post_appen(t,'4_cMul Elements');
write_cmat(f,c);
tic; c = bsxfun(@times, ddd,iii );
total = total + toc; f = post_appen(t,'5_cMul Elements');
write_cmat(f,c);
% 1 4
tic; c = bsxfun(@times, idi,did );
total = total + toc; f = post_appen(t,'6_cMul Elements');
write_cmat(f,c);
tic; c = bsxfun(@times, dii,idd );
total = total + toc; f = post_appen(t,'7_cMul Elements');
write_cmat(f,c);
tic; c = bsxfun(@times, idi, ddd);
total = total + toc; f = post_appen(t,'8_cMul Elements');
write_cmat(f,c);
% 1 1
tic; c = bsxfun(@rdivide,dii ,ddi );
total = total + toc; f = post_appen(t,'9_cDiv Elements');
write_cmat(f,c);
tic; c = bsxfun(@rdivide, ddi,idi );
total = total + toc; f = post_appen(t,'10_cDiv Elements');
write_cmat(f,c);
tic; c = bsxfun(@rdivide, ddi,dii );
total = total + toc; f = post_appen(t,'11_cDiv Elements');
write_cmat(f,c);
tic; c = bsxfun(@rdivide, ddi,ddi );
total = total + toc; f = post_appen(t,'12_cDiv Elements');
write_cmat(f,c);

y = total;

end

