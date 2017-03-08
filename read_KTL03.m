function D=read_KTL03(filename)

GT={'1','2','3'};
LR={'L','R'};

out_struct.geolocation=struct(...
    'segment_dist_x',[],...
    'segment_id', [], ...
    'surf_type', []);
out_struct.heights=struct(...
    'delta_time', [], ...
    'dist_ph', [], ...
    'dist_ph_across', [], ...
    'dist_ph_along', [], ...
    'h_ph', [], ...
    'lat_ph', [], ...
    'lon_ph', [], ...
    'pce_mframe_cnt', [], ...
    'signal_conf_ph', []);

D=repmat(out_struct, [6,1]);

%fprintf(1,'%s\n beam\t N_ph \t N_signal_conf_ph\n', filename);
for kT=1:length(GT)
    for kB=1:length(LR)
        beam=(kT-1)*2+kB;
        GT_grp=sprintf('/GT%s%s', GT{kT}, LR{kB});
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


for beam=1:length(D)
    % extract the ice-sheet column for the signal confidence
    D(beam).heights.signal_conf_ph=D(beam).heights.signal_conf_ph(4,:)';
    % extract the segment distance values using pce_mframe_cnt
    D(beam).heights.dist_seg=D(beam).heights.dist_ph(D(beam).heights.pce_mframe_cnt);
    D(beam).heights=rmfield(D(beam).heights,'dist_ph');
end

% example plot:

if false
    thefile='/Volumes/ice1/ben/sdt/KTL03/KTL03s_for_ben/ATL03_sim_photons_allbeams_seg2_rgt696_v5.h5';
    D=read_KTL03(thefile);
    for kB=1:6
        temp=D(kB).heights;
        temp=index_struct(temp, 1:100:length(temp.h_ph));
        figure(kB); clf; hold on; 
        for val=0:4
            els=temp.signal_conf_ph==val;
            plot(temp.dist_seg(els)+temp.dist_ph_along(els), temp.h_ph(els),'.');
        end
    end
    
end

