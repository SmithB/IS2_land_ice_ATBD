function h=ATL11_dumbell_plot(y, h, rep, beam, varargin)

[~, ind]=sort(rep+y/1e6);
ind=reshape(ind, 2, length(y(:))/2);
if nargin > 3
    plot(y(ind), h(ind), varargin{:});
else
    plot(y_ind, h(ind),'linestyle','-');
end

ht=text(y(:), h(:), num2str(rep(:)));

