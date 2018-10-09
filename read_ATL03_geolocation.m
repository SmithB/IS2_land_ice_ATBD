function [D, dist_for_segment]=read_ATL03_geolocation(filename, pairs, params, geoloc_full)

 
if ~exist('pairs','var')
    pairs=1:3;
end

beams=sort([2*pairs(:)-1; 2*pairs(:)]);

if ~isfield(params,'index_range')
    for kP=1:length(pairs)
        for kB=1:2
            params.index_range{kB, pairs(kP)}=[1 Inf];
        end
    end
end
dist_for_segment=sparse([]);
 
GT={'1','2','3'};
LR={'l','r'};
if isfield(params,'seg_info_only') && params.seg_info_only
    out_struct=struct(...
        'segment_dist_x',[],...
        'segment_id', [], ...
        'ph_index_beg', [], ...
        'segment_ph_cnt', []);
else
    out_struct=struct(...
    'segment_dist_x',[],...
    'segment_id', [], ...
    'ph_index_beg', [], ...
    'delta_time', [], ...
    'segment_ph_cnt', [], ...
    'sigma_across', [], ...
    'sigma_along', []);
end
% write the empty struct to D (in case we have to return early), and
% add
D=repmat(setfield(out_struct, 'geoloc_rec_index', []), [max(beams),1]);


if isfield(params,'seg_range')
    % a range of segments were requested.  We need to have already read the
    % geoloc file to do this-- if the geoloc_full variable wasn't
    % specified, read teh file
    if ~exist('geoloc_full','var')
        geoloc_full=read_ATL03_geolocation(filename, pairs, struct('seg_info_only',true));
    end
    found_indices=false;
    for kP=pairs
        for kB=1:2
            % find the range of indices we need to read to get requested
            % seg ids
            beam=(kP-1)*2+kB;
            i0=find(geoloc_full(beam).segment_id >= params.seg_range(1), 1, 'first');            
            i1=find(geoloc_full(beam).segment_id <= params.seg_range(2), 1, 'last');
            if isempty(i0) || isempty(i1) || i0>i1
                % Oops.  Didn't find signals.  Return empty arrays
                continue
            end
            params.index_range{kB, kP}=geoloc_full(beam).geoloc_rec_index([i0 i1]);
            found_indices=true;
        end
    end
    if ~found_indices; return; end
end

 
for kP=pairs(:)' 
    for kB=1:length(LR)
        % read in the parameters for each beam
        beam=(kP-1)*2+kB;
        GT_grp=sprintf('/gt%s%s', GT{kP}, LR{kB});
        F1=fieldnames(out_struct);
        for k1=1:length(F1)
            fieldName=[GT_grp,'/geolocation/', F1{k1}];           
            D(beam).(F1{k1})=h5read(filename, fieldName, params.index_range{kB, kP}(1), diff(params.index_range{kB, kP})+1);
        end
        D(beam).geoloc_rec_index=(1:length(D(beam).(F1{k1})))';
        D(beam)=index_struct(D(beam), D(beam).ph_index_beg ~=0);
        
        if isfield(D,'segment_dist_x') && nargout>1
            dist_for_segment{beam}=sparse(double(D(beam).segment_id), ones(size(D(beam).segment_id)), D(beam).segment_dist_x);
        end
    end
end
 