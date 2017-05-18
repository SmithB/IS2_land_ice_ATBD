
N_pulses=57*10;
BGR_vals=logspace(5.5, 7, 20);
Htot=200;
N_chan=[4, 16];  %weak beam (but doesn't matter much for noise)
load WF_est;   

 for kBG=1:1%length(BGR_vals);
    this_file=sprintf('noise_test_files/BGR=%3.2f.h5', BGR_vals(kBG)/1e6);
    
    clear D2a;
    for kB=1:2
        [D2, params(kB)]=make_ATL03_data(N_pulses, N_chan(kB), 0, 0, WF, BGR_vals(kBG), Htot, 0);
        D2=index_struct(D2, D2.detected);
        % delete problem parameters
        non_atl03_parameters={'x0','zground','xground','SigNoise' };
        for k3=1:length(non_atl03_parameters)
            D2=rmfield(D2, non_atl03_parameters{k3});
        end
        
        D2.segment_number=floor(D2.pulse_num/29);
        D2.x_RGT=D2.pulse_num*0.7;
        D2.y_RGT=90*(kB-0.5)*ones(size(D2.h));
        D2.time=D2.pulse_num/1e4 + D2.t_ph;
        D2.ph_class=zeros(size(D2.h));
        [D2.beam, D2.track]=deal(ones(size(D2.h))*kB);
        
        
        D2a(kB)=D2;
     end
    if exist(this_file,'file');;
        delete(this_file)
    end
    make_test_data('write_test_file', this_file, D2a, params);
end