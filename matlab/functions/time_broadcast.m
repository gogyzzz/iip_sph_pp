function y = time_broadcast( iii,iid,idi,idd,dii,did,ddi,ddd )

tic;
c = bsxfun(@plus, iid, idd );
c = bsxfun(@plus, iid,did );
c = bsxfun(@plus, idd,iid );
% 4 1
c = bsxfun(@plus, did,iii );
c = bsxfun(@times, iid, ddi);
c = bsxfun(@times, ddd,iii );
% 1 4
c = bsxfun(@times, idi,did );
c = bsxfun(@times, dii,idd );
c = bsxfun(@times, idi, ddd);
% 1 1
c = bsxfun(@rdivide,dii ,ddi );
c = bsxfun(@rdivide, ddi,idi );
c = bsxfun(@rdivide, ddi,dii );
c = bsxfun(@rdivide, ddi,ddi );
y= toc;

end


