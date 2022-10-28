function [D_06, residual_histogram]=read_ASAS_ATL06(thefile, pairs)

if ~exist('pairs','var')
    pairs=[1 2 3];
end

beams={'gt1l','gt1r','gt2l','gt2r','gt3l','gt3r'};
% altimetry:
fields.None={'segment_id','h_li','h_li_sigma', 'delta_time','atl06_quality_summary' ,'latitude','longitude'};
fields.ground_track={'x_atc','y_atc'};
fields.fit_statistics={'dh_fit_dx', 'dh_fit_dy', 'dh_fit_dx_sigma', 'n_fit_photons','n_seg_pulses', ...
    'h_robust_sprd', 'w_surface_window_final','signal_selection_source','snr_significance', 'h_mean'};

fields.bias_correction={'fpb_n_corr','tx_med_corr','fpb_med_corr','fpb_mean_corr'};
fields.geophysical={'bckgrd','dac','tide_ocean'};
for kP=pairs
    for kB=1:2
        for group=fieldnames(fields)'
             if strcmp(group{1},'None')
                 thegroup=sprintf('/%s/land_ice_segments', beams{2*(kP-1)+kB});
             else
                 thegroup=sprintf('/%s/land_ice_segments/%s', beams{2*(kP-1)+kB}, group{1} );
             end
             for thefield=fields.(group{1})
                thevar=sprintf('%s/%s', thegroup, thefield{1});
                D_06(kP).(thefield{1})(:, kB)=read_h5_var(thefile, thevar);
             end           
        end
 
        if nargout > 1
            thegroup=sprintf('/%s/residual_histogram', beams{2*(kP-1)+kB});
            try
                II=h5info(thefile, thegroup);
            catch
                continue
            end
            for kF=1:length(II.Datasets);
                residual_histogram(kB, kP).(II.Datasets(kF).Name)=read_h5_var(thefile, [thegroup,'/', II.Datasets(kF).Name]);
            end
        end
    end
end


for kP=pairs;
    D_06(kP).h_li(D_06(kP).h_li>1e4)=NaN;
    D_06(kP).n_seg_pulses(D_06(kP).n_seg_pulses==0)=NaN;
end

%--------------------------------
function x=read_h5_var(thefile, thevar);

x=double(h5read(thefile, thevar));
try
    temp=double(h5readatt(thefile,thevar,'_FillValue'));
    x(x==temp)=NaN;
end

