% unit test 1.

% input: 
% h_ph  photon heights
% segment_number : the segment number we want
% ph_class : photon classification
%
%  For visualization:
% x_RGT : along-track coordinate (for visualization)




% pick a data set
test_data_dir=['/home/ben/Dropbox/projects/IS2_ATBD/Combined_corr_test_data_May_4_2017/'];


% pick data files
D2_file='D2/Rough=4.80e-01_Rsurf=1.00e+00_BGR=1.00e+07.h5';
D3_file='D3/Rough=4.80e-01_Rsurf=1.00e+00_BGR=1.00e+07.h5';

% select segment number 7:
m=7;
x0=20*(m-1);
params.sigma_x=5;
params.sigma_pulse=0.68e-9;

%h5disp([test_data_dir, test_file],'/','min')

D2_group='/photon/weak';
D2_fields={'detected','h','ph_class','segment_number', 'x_RGT', 'BGR', 'pulse_num'};
% read the input data:
for kf=1:length(D2_fields)
    D2.(D2_fields{kf})=h5read([test_data_dir, D2_file], [D2_group,'/', D2_fields{kf}]);
end
% read the output data:
D3_fields={'seg_count','segment_id','m0','N_initial','w_surface_window_initial','x_RGT'};
for kf=1:length(D3_fields)
    D3.(D3_fields{kf})=h5read([test_data_dir, D3_file], ['/', D3_fields{kf}]);
end

% select ATL03 segments m-1 and m
AT_els=[D2.segment_number==m-1 | D2.segment_number==m];
D2sub=index_struct(D2, AT_els & D2.detected);

% call choose_ground_strategy
[initial_fit_els, signal_selection_source, signal_selection_status_confident, signal_selection_status_all]=...
    ATLAS_L3a_proc_ATBD('choose_ground_strategy', D2sub);

if signal_selection_source < 1  % enough confident PE found to determine slope and height of the bin
    initial_Hwin_min=3;
else %
    initial_Hwin_min=10;
end


params.pulses_per_seg=length(unique(D2sub.pulse_num));
% Determine the window size that corresponds to the selected PE
[initial_fit_els, w_surface_window_initial, h_initial, ~]=...
    ATLAS_L3a_proc_ATBD('initial_at_fit', D2sub.x_RGT, D2sub.h, initial_fit_els, x0, median(D2sub.BGR), initial_Hwin_min, params);
        
D3_row=find(D3.segment_id(:,1)==m);
D3_col=1;
D3.N_initial(D3_row, D3_col)


        
% 








