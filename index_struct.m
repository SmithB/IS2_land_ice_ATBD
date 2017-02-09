function S1=index_struct(S, els, f);

% S1=index_struct(S, els, f);
%
% selects the elements (or rows) in the fields in structure S corresponding
% to index array 'els'.  if a cell array of field names is provided in 'f', 
% the output structure contains only those fields from S.


if ~exist('f','var');
    f=fieldnames(S);
end
for k=1:length(f);
    eval(['F=size(S.', f{k},');']);
    if any(F==0)
        eval(['S1.', f{k},'=[];']);
        continue
    end
    if F(2)>1 && length(els(:))==F(1);
        eval(['S1.', f{k},'=S.',f{k},'(els,:);']);
    else
        eval(['S1.', f{k},'=S.', f{k},'(els);']);
    end      
end
