function D=read_KTL03(filename)

GT={'1','2','3'};
LR={'L','R'};

out_struct=struct('delta_time', [], 'dist_ph_across',[], 'dist_ph_along',[], 'lat',[], 'lon',[],'segment_dist_x',[]);
 
D=repmat(out_struct, [6,1]);

%fprintf(1,'%s\n beam\t N_ph \t N_signal_conf_ph\n', filename);
for kT=1:length(GT)
    for kB=1:length(LR)
        beam=(kT-1)*2+kB;
        GT_grp=sprintf('/GT%s%s', GT{kT}, LR{kB});
        F0=fieldnames(out_struct);
        for k0=1:length(F0)           
            fieldName=[GT_grp,'/',F0{k0}];
            D(beam).(F0{k0}) =h5read(filename, fieldName);          
        end
        %fprintf(1, '%s:\t%8.0d\t%8.0d\n',  GT_grp, length(D(beam).heights.h_ph), size(D(beam).heights.signal_conf_ph,2));
    end
end
