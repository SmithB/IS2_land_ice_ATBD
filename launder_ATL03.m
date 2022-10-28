function D2a=launder_ATL03(D2a)

delete_fields={};
if ~isfield(D2a,'segment_id')
    if isfield(D2a,'seg_num')
        for k=1:length(D2a)
            D2a(k).segment_id=D2a(k).seg_num;
        end
        delete_fields{end+1}='seg_num';
    elseif isfield(D2a,'segment_number')
        for k=1:length(D2a)
            D2a(k).segment_id=D2a(k).segment_number;
        end
        delete_fields{end+1}='segment_number';
    end
end

if ~isfield(D2a,'h_ph')
    if isfield(D2a,'h')
        for k=1:length(D2a)
            D2a(k).h_ph=D2a(k).h;
        end
        delete_fields{end+1}='h';
    elseif isfield(D2a, 'elev')
        for k=1:length(D2a)
            D2a(k).h_ph=D2a(k).elev;
        end
        delete_fields{end+1}='elev';
    end
end

if ~isfield(D2a,'signal_conf_ph')
    for k=1:length(D2a)
        D2a(k).signal_conf_ph=D2a(k).ph_class;
    end
    delete_fields{end+1}='ph_class';
end

for k=1:length(delete_fields)
    D2a=rmfield(D2a, delete_fields{k});
end
        
        
        