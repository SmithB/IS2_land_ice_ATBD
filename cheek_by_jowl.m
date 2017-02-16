function AX=cheek_by_jowl(rows, cols, box, delta);

%AX=cheek_by_jowl(rows, cols, box, delta);
% generates a set of (rows, cols) axes that
% nearly fill a rectangle defined by 'box'.
% 
% box is in normalized coordinates in the current
% figure window, in the standard matlab format
% [min_x min_y delta_x delta_y]
%
% the optional parameter 'delta' specifies
% the spacing between x and y axes, respectively, in normalized units.
% the default value is [0.01 0.01] 

if ~exist('delta','var');
    delta=[0.01 0.01];
end

if length(delta)==1
    delta=[delta delta];
end

dx=box(3)/cols; dy=box(4)/rows; 

if dx==box(3)
    x0=box(1);
else
    x0=box(1):dx:(box(1)+box(3)-dx);
end
if dy==box(4);
    y0=box(2);
else
    y0=box(2):dy:(box(2)+box(4)-dy);
end
y0=flipud(y0(:));

for row=1: rows; 
    for col=1:cols; 
        AX(row, col)=axes('position', [x0(col)+delta(1), y0(row)+delta(2), dx-delta(1)/2, dy-delta(2)/2]);
    end
end

set(AX(1:end-1,:),'xtick', []); set(AX(:,2:end-1), 'ytick', []);
set(AX(:,end),'yaxislocation', 'right');



