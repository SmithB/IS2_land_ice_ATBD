function B=bilinear_basis_functions(x, y, x0, y0)

x=double(x); y=double(y);

[xg, x0g]=ndgrid(x(:), x0);
[yg, y0g]=ndgrid(y(:), y0);

if length(x0) > 1
    Bx=sparse([], [], [], size(xg, 1), size(xg,2), min([length(x0) 4])*length(x));
    delta_norm=(xg(:, 1:end-1)-x0g(:, 1:end-1))./(x0g(:, 2:end)-x0g(:, 1:end-1));
    Bx(:, 2:end)=delta_norm.*(delta_norm >0 & delta_norm < 1);
    Bx(x0g==x0(end) & xg==x0(end))=1;
    delta_norm=(x0g(:, 2:end)-xg(:, 1:end-1))./(x0g(:, 2:end)-x0g(:, 1:end-1));
    Bx(:, 1:end-1)=Bx(:, 1:end-1)+ delta_norm.*(delta_norm >0 & delta_norm <= 1);
else
    Bx=ones(size(x(:))); 
end

if length(y0) > 1
    By=sparse([], [],[], size(yg, 1), size(yg,2), min([length(y0),4])*length(x));
    delta_norm=(yg(:, 1:end-1)-y0g(:, 1:end-1))./(y0g(:, 2:end)-y0g(:, 1:end-1));
    By(:, 2:end)=delta_norm.*(delta_norm >0 & delta_norm < 1);
    delta_norm=(y0g(:, 2:end)-yg(:, 1:end-1))./(y0g(:, 2:end)-y0g(:, 1:end-1));
    By(:, 1:end-1)=By(:, 1:end-1)+ delta_norm.*(delta_norm >0 & delta_norm <= 1);
    By(y0g==y0(end) & yg==y0(end))=1;
else
    By=ones(size(x(:)));
end

% count=0;
% B=sparse(length(x(:)), length(x0)*length(y0));
% for kX=1:length(x0);
%     if any(Bx(:, kX))
%         for kY=1:length(y0);
%             count=count+1;
%             if any(By(:, kY));
%                 B(:, count)=Bx(:,kX).*By(:, kY);
%             end
%         end
%     else
%         count=count+length(y0);
%     end
% end

[Byc, Bxc]=ind2sub([length(y0), length(x0)], 1:(length(x0)*length(y0)));
% want this: B=Bx(:, Bxc).*By(:, Byc);
% sad to say, it takes too much memory. 
% so loop over the coumns of Bx and By like this:
col_ind=[1:1000:length(Bxc) length(Bxc)+1];
B=sparse(length(x(:)), length(x0)*length(y0));
for k=1:length(col_ind)-1;
    these_cols=col_ind(k):col_ind(k+1)-1;
    B(:, these_cols)= Bx(:, Bxc(these_cols)).*By(:, Byc(these_cols));
end




