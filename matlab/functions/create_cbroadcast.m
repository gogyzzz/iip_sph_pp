function create_cbroadcast( a0,a1,a2,a3,a4,a5,a6,a7 )

iii = read_cmat(a0);
iid = read_cmat(a1);
idi = read_cmat(a2);
idd= read_cmat(a3);
dii = read_cmat(a4);
did= read_cmat(a5);
ddi= read_cmat(a6);
ddd = read_cmat(a7);

% 4 4
c = bsxfun(@plus, iid, idd );
write_mat('data/cBroadcast_0_add.bin',c);
c = bsxfun(@plus, iid,did );
write_mat('data/cBroadcast_1_add.bin',c);
c = bsxfun(@plus, idd,iid );
write_mat('data/cBroadcast_2_add.bin',c);
% 4 1
c = bsxfun(@plus, did,iii );
write_mat('data/cBroadcast_3_add.bin',c);
c = bsxfun(@times, iid, ddi);
write_mat('data/cBroadcast_4_mul.bin',c);
c = bsxfun(@times, ddd,iii );
write_mat('data/cBroadcast_5_mul.bin',c);
% 1 4
c = bsxfun(@times, idi,did );
write_mat('data/cBroadcast_6_mul.bin',c);
c = bsxfun(@times, dii,idd );
write_mat('data/cBroadcast_7_mul.bin',c);
c = bsxfun(@times, idi, ddd);
write_mat('data/cBroadcast_8_mul.bin',c);
% 1 1
c = bsxfun(@rdivide,dii ,ddi );
write_mat('data/cBroadcast_9_div.bin',c);
c = bsxfun(@rdivide, ddi,idi );
write_mat('data/cBroadcast_10_div.bin',c);
c = bsxfun(@rdivide, ddi,dii );
write_mat('data/cBroadcast_11_div.bin',c);
c = bsxfun(@rdivide, ddi,ddi );
write_mat('data/cBroadcast_12_div.bin',c);



end

