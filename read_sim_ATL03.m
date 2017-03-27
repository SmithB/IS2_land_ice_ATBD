function [D, params]=read_sim_ATL03(filename, pairs)

if ~exist('pairs','var')
    pairs=1:3;
end

beams=sort([2*pairs(:)-1; 2*pairs(:)]);


GT={'1','2','3'};
LR={'l','r'};

out_struct.geolocation=struct(...
    'segment_dist_x',[],...
    'segment_id', [], ...
    'surf_type', [], ...
    'ph_index_beg', [], ...
    'delta_time', [], ...
    'segment_ph_cnt', []);
out_struct.heights=struct(...
    'delta_time', [], ...
    'dist_ph_across', [], ...
    'dist_ph_along', [], ...
    'h_ph', [], ...
    'lat_ph', [], ...
    'lon_ph', [], ...
    'pce_mframe_cnt', [], ...
    'ph_id_count', [], ...
    'signal_conf_ph', []);
out_struct.bckgrd_atlas=struct('bckgrd_rate', [],'pce_mframe_cnt',[]);

D=repmat(out_struct, [max(beams),1]);

%fprintf(1,'%s\n beam\t N_ph \t N_signal_conf_ph\n', filename);
for kT=pairs(:)'%length(GT)
    for kB=1:length(LR)
        beam=(kT-1)*2+kB;
        GT_grp=sprintf('/gt%s%s', GT{kT}, LR{kB});
        F0=fieldnames(out_struct);
        for k0=1:length(F0)
            F1=fieldnames(out_struct.(F0{k0}));
            for k1=1:length(F1)
                fieldName=[GT_grp,'/',F0{k0},'/', F1{k1}]; 
                D(beam).(F0{k0}).(F1{k1})=h5read(filename, fieldName);
            end
        end  
        %fprintf(1, '%s:\t%8.0d\t%8.0d\n',  GT_grp, length(D(beam).heights.h_ph), size(D(beam).heights.signal_conf_ph,2));
    end
end


for beam=beams(:)' %length(D)
    % extract the ice-sheet column for the signal confidence
    D(beam).heights.signal_conf_ph=D(beam).heights.signal_conf_ph(4,:)';
    
    [D(beam).heights.ph_seg_num, ...
        D(beam).heights.x_RGT]=deal(NaN(size(D(1).heights.h_ph)));
    for k=1:length(D(beam).geolocation.ph_index_beg)-1
        ind_range=[D(beam).geolocation.ph_index_beg(k), D(beam).geolocation.ph_index_beg(k+1)-1];
        if all(ind_range~=0)
            D(beam).heights.ph_seg_num(ind_range(1):ind_range(2))=k;
        end
    end
    D(beam).heights.ph_seg_num(D(beam).geolocation.ph_index_beg(end):end)=k+1;
    these=isfinite(D(beam).heights.ph_seg_num); 
    D(beam).heights.x_RGT(these)=D(beam).geolocation.segment_dist_x(D(beam).heights.ph_seg_num(these));%+D(beam).heights.dist_ph_along(these);
    good=isfinite(D(beam).heights.pce_mframe_cnt & D(beam).heights.pce_mframe_cnt & D(beam).heights.pce_mframe_cnt < length(D(beam).bckgrd_atlas.bckgrd_rate));
    D(beam).BGR(good)=D(beam).bckgrd_atlas.bckgrd_rate(D(beam).heights.pce_mframe_cnt(good));
end


if nargout>1
    for k=1:length(D)
        params(k).WF.t=h5read(filename,'/atlas_impulse_response/tep/beam_3/histogram/tep_hist_time');
        params(k).WF.p=h5read(filename,'/atlas_impulse_response/tep/beam_3/histogram/tep_hist');
    end
end

% example plot:

if false
    thefile='/Volumes/ice1/ben/sdt/KTL03/KTL03s_for_ben/ATL03_sim_photons_allbeams_seg2_rgt696_v5.h5';
    D=read_sim_ATL03(thefile);
     for kB=1:2
        temp=D(kB).heights;
        temp=index_struct(temp, 1:100:length(temp.h_ph));
        figure(kB); clf; hold on; 
        for val=0:4
            els=temp.signal_conf_ph==val;
            plot(temp.x_RGT(els), temp.h_ph(els),'.');
        end
     end
     
end

