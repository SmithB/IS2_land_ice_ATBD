function [H, D]=read_ATL03_photon_data(filename, pairs, params, geoloc)

beams=sort([2*pairs(:)-1; 2*pairs(:)]);

GT={'1','2','3'};
LR={'l','r'};

if ~exist('params','var'); params=struct(); end

if ~isfield(params, 'fields')
    out_struct=struct(...
        'delta_time', [], ...
        'dist_ph_across', [], ...
        'dist_ph_along', [], ...
        'h_ph', [], ...
        'lat_ph', [], ...
        'lon_ph', [], ...
        'pce_mframe_cnt', [], ...
        'ph_id_count', [], ...
        'ph_id_pulse', [], ...
        'ph_id_channel',[], ...
        'signal_conf_ph', []);
end
bckgrd_struct=struct('bckgrd_rate', [],'pce_mframe_cnt',[],'delta_time', []);

if ~isfield(params,'signal_conf_type')
    params.signal_conf_type=4;
end
 
if ~exist('geoloc','var')
    geoloc=read_ATL03_geolocation(filename, pairs, params);
end

% make the photon count range (tells which photons to read
for kP=pairs(:)'
    for kB=1:length(LR)
        beam=(kP-1)*2+kB;
        if ~isempty(geoloc(beam).ph_index_beg)
            ph_cnt_range{kB, kP}=[geoloc(beam).ph_index_beg(1) geoloc(beam).ph_index_beg(end)+int64(geoloc(beam).segment_ph_cnt(end))];
        else
            ph_cnt_range{kB, kP}=[];
        end
    end
end

   
% read the photon data beam by beam
for kT=pairs(:)'%length(GT)
    for kB=1:length(LR)
        beam=(kT-1)*2+kB;
        GT_grp=sprintf('/gt%s%s', GT{kT}, LR{kB});
        
        F1=fieldnames(out_struct);
        for k1=1:length(F1)
            fieldName=[GT_grp,'/heights/', F1{k1}];
            if isempty(ph_cnt_range{kB, kT}); continue; end
            start= double(ph_cnt_range{kB, kT}(1));
            count= double(diff(ph_cnt_range{kB, kT}));
            if ~strcmp(F1{k1},'signal_conf_ph')
                D(beam).(F1{k1})=h5read(filename, fieldName, start, count);
            else
                % for signal_conf, need to read just one row out of the
                % dataset
                D(beam).(F1{k1})=h5read(filename, fieldName, [params.signal_conf_type, start], [1,count])';
            end
        end
        try
            %... and read the background info
            F1=fieldnames(bckgrd_struct);
            bckgrd(beam)=bckgrd_struct;
            for kF=1:length(F1)
                bckgrd(beam).(F1{kF})=h5read(filename,[GT_grp,'/bckgrd_atlas/' F1{kF}]);
            end
        catch
            fprintf(1,'background info missing for %s', GT_grp);
        end
    end
end

% map the geoloc parameters to the photon parametes
for kT=pairs(:)'%length(GT)
    for kB=1:length(LR)
        beam=(kT-1)*2+kB;
        if isempty(ph_cnt_range{kB, kT}); continue; end

        
        % extract the ice-sheet column for the signal confidence
        H(beam).signal_conf_ph=D(beam).signal_conf_ph;
        
        [H(beam).ph_seg_num, ...
            H(beam).x_RGT, ...
            H(beam).seg_dist_x, ...
            H(beam).BGR, ...
            H(beam).sigma_across, ...
            H(beam).sigma_along, ...
            H(beam).segment_id]=deal(NaN(size(D(beam).h_ph)));
        % we have read photons starting with ph_cnt_range(beam)
        offset=ph_cnt_range{kB, kT}(1)-1;
        for k=1:length(geoloc(beam).ph_index_beg)
            ind_range=geoloc(beam).ph_index_beg(k)+int64([0 geoloc(beam).segment_ph_cnt(k)-1]);
            if all(ind_range>0)
                H(beam).ph_seg_num((ind_range(1):ind_range(2))-offset)=k;
            end
        end
        
        these=isfinite(H(beam).ph_seg_num) & H(beam).ph_seg_num > 0;
        H(beam).seg_dist_x(these)=geoloc(beam).segment_dist_x(H(beam).ph_seg_num(these));
        H(beam).sigma_across(these)=geoloc(beam).sigma_across(H(beam).ph_seg_num(these));
        H(beam).sigma_along(these)=geoloc(beam).sigma_along(H(beam).ph_seg_num(these));
        
        
        H(beam).x_RGT(these)=geoloc(beam).segment_dist_x(H(beam).ph_seg_num(these))+double(D(beam).dist_ph_along(these));
        H(beam).dist_ph_across=D(beam).dist_ph_across;
        
        good=isfinite(H(beam).x_RGT) & D(beam).pce_mframe_cnt>0; %& H(beam).pce_mframe_cnt < length(D(beam).bckgrd_atlas.bckgrd_rate);
        try
            H(beam).BGR(good)=interp1(bckgrd(beam).delta_time, bckgrd(beam).bckgrd_rate, D(beam).delta_time(good));
        catch
            H(beam).BGR(good)=NaN(size(D(beam).delta_time(good)));
        end
        H(beam).pulse_num=200*double(D(beam).pce_mframe_cnt)+double(D(beam).ph_id_pulse);
        H(beam).beam=zeros(size(D(beam).h_ph))+beam;
        H(beam).h_ph=double(D(beam).h_ph);
        H(beam).delta_time=D(beam).delta_time;
        H(beam).lat_ph=D(beam).lat_ph;
        H(beam).lon_ph=D(beam).lon_ph;
        
        H(beam).segment_id(these)=geoloc(beam).segment_id(H(beam).ph_seg_num(these));
        H(beam).ph_seg_num=H(beam).ph_seg_num+double(offset);
        H(beam).ph_id_channel=D(beam).ph_id_channel;
    end
end

