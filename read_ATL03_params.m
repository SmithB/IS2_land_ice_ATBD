function params=read_ATL03_params(filename, beams)

if ~exist('pairs','var')
    pairs=1:3;
end
 
GT={'1','2','3'};
LR={'l','r'};     

% read which TEP to use with which spot (laser beam?? check this)
try
    TEP_spots=h5read(filename,'/ancillary_data/tep/tep_valid_spot');
catch
    TEP_spots=ones(1,6);
end
% read the TEP
TEP_groups={'pce1_spot1', 'pce2_spot3'};
for kG=1:length(TEP_groups)
    WF(kG).t=h5read(filename, ['/atlas_impulse_response/', TEP_groups{kG},'/tep_histogram/tep_hist_time']);
    WF(kG).p=h5read(filename, ['/atlas_impulse_response/', TEP_groups{kG},'/tep_histogram/tep_hist']);
    WF(kG)=proc_TEP([], [], WF(kG));
end

for kT=pairs(:)'%length(GT)
    for kB=1:length(LR)
        beam=(kT-1)*2+kB;;
        GT_grp=sprintf('/gt%s%s', GT{kT}, LR{kB});      

         if strcmp(deblank(h5readatt(filename, GT_grp,'atlas_beam_type')),'strong')
            params(beam).N_det=16;
        else
            params(beam).N_det=4;
        end
        params(beam).RGT=h5read(filename,'/ancillary_data/start_rgt');
        params(beam).orbit=h5read(filename,'/ancillary_data/start_orbit');
        params(beam).cycle=h5read(filename,'/ancillary_data/start_cycle');
        params(beam).GT=beam;
        params(beam).spot_number=str2double(deblank(h5readatt(filename, GT_grp,'atlas_spot_number')));
        params(beam).PT=str2double(kT);     
        params(beam).WF=WF(TEP_spots(params(beam).spot_number));
         
    end
end
 


