function create_broadcast( a0,a1,a2,a3,a4,a5,a6,a7 )

iii = read_mat(a0);
iid = read_mat(a1);
idi = read_mat(a2);
idd= read_mat(a3);
dii = read_mat(a4);
did= read_mat(a5);
ddi= read_mat(a6);
ddd = read_mat(a7);

% 4 4
c = bsxfun(@plus, iid, idd );
write_mat('data/Broadcast_0_add.bin',c);
c = bsxfun(@plus, iid,did );
write_mat('data/Broadcast_1_add.bin',c);
c = bsxfun(@plus, idd,iid );
write_mat('data/Broadcast_2_add.bin',c);
% 4 1
c = bsxfun(@plus, did,iii );
write_mat('data/Broadcast_3_add.bin',c);
c = bsxfun(@times, iid, ddi);
write_mat('data/Broadcast_4_mul.bin',c);
c = bsxfun(@times, ddd,iii );
write_mat('data/Broadcast_5_mul.bin',c);
% 1 4
c = bsxfun(@times, idi,did );
write_mat('data/Broadcast_6_mul.bin',c);
c = bsxfun(@times, dii,idd );
write_mat('data/Broadcast_7_mul.bin',c);
c = bsxfun(@times, idi, ddd);
write_mat('data/Broadcast_8_mul.bin',c);
% 1 1
c = bsxfun(@rdivide,dii ,ddi );
write_mat('data/Broadcast_9_div.bin',c);
c = bsxfun(@rdivide, ddi,idi );
write_mat('data/Broadcast_10_div.bin',c);
c = bsxfun(@rdivide, ddi,dii );
write_mat('data/Broadcast_11_div.bin',c);
c = bsxfun(@rdivide, ddi,ddi );
write_mat('data/Broadcast_12_div.bin',c);



end

