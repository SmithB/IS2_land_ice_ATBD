function [dM, dCtr, syn_wf]=correct_for_TX_shape( t_TX, TX, HW, sigma_rx, SNR)

% inputs:
% t_tx : Time scale for the transmit pulse
% TX: power for the transmit pulse
% HW: height of the window, in same units as t_tx
% sigma_rx: width of the rx pulse (same units as t_tx)
% SNR: signal-to-noise rate within the RX window

% returns: dM, dCtr: height offsets of the median and mean of the rx pulse
% WRT the centroid of the TX pulse.  Add these to the uncorrected mean and
% median heights

sigma_tx=(wf_percentile(t_TX, TX, 0.84)-wf_percentile(t_TX, TX, 0.16))/2;
dt=t_TX(2)-t_TX(1);

% generate a WF that looks like the tx pulse convolved with the pulse-spreading function
% broaden the TX waveform if it's narrower than the RX waveform
if sigma_rx > sigma_tx
    sigma_G=sqrt(sigma_rx.^2-sigma_tx.^2);
    Gsamps=ceil(sigma_G/dt);
    G=gaussian(dt*(-4*Gsamps:4*Gsamps), 0, sigma_G);
    t1=[(-4*Gsamps:-1)*dt+t_TX(1) t_TX(:)' (1:4*Gsamps)*dt+t_TX(end)]; t1=t1(:);
else
    G=1;
    t1=t_TX(:);
end
TX_WF_broadened=conv(TX(:), G(:));
 
% find the 16th and 84th percentiles of the broadened WF, use these to truncate it
LCH=wf_percentile(t1, TX_WF_broadened, [0.16 0.5 0.84]);

if exist('SNR', 'var')
    N_add_early=ceil((1.25*HW/2-(LCH(2)-t1(1)))/dt);
    if N_add_early > 0
        % add samples to the front of the WF
        t1=[t1(1)+(-N_add_early:-1)'*dt ; t1];
        TX_WF_broadened=[zeros(N_add_early,1); TX_WF_broadened];
    end
    N_add_late=ceil((1.25*HW/2-(t1(end)-LCH(2)))/dt);
    if N_add_late > 0
        % add samples to the front of the WF
        t1=[t1; t1(end)+(1:N_add_late)'*dt ];
        TX_WF_broadened=[TX_WF_broadened; zeros(N_add_late,1)  ];
    end
       
    TX_WF_broadened=TX_WF_broadened/sum(TX_WF_broadened)+1/SNR/(HW/dt);

    % iterate to find the centroid
    last_t_ctr=t1(1);
    % begin centered on the median
    t_ctr=LCH(2);
    while abs(t_ctr-last_t_ctr) > 1e-4/1.5e8
        els=abs(t1-t_ctr)<HW/2;
        last_t_ctr=t_ctr;
        t_ctr=sum(TX_WF_broadened(els).*t1(els))/sum(TX_WF_broadened(els));
    end
    TXBm0=t_ctr;
else
    % truncate around the median
    TXB_els=abs(t1-LCH(2)) < HW/2;%1.5*(LCH(3)-LCH(1));
    
    %find the mean of the RX waveform (edited by twice the window)
    %TXBm0=wf_percentile(t1(TXB_els), TX_WF_broadened(TXB_els), 0.5);
    TXBm0=sum(t1(TXB_els).*TX_WF_broadened(TXB_els))./sum(TX_WF_broadened(TXB_els));
    els=TXB_els;
end

%truncate the synthetic WF around its edited mean by HW
TXB_els_HW=abs(t1-TXBm0) < HW/2;
TXBm=wf_percentile(t1(TXB_els_HW), TX_WF_broadened(TXB_els_HW), 0.5);
 
%dM=TXm-RXm;
TimeToH=1.5e8;
dM=TXBm*TimeToH;
dCtr=sum(t1(TXB_els_HW).*TX_WF_broadened(TXB_els_HW))./sum(TX_WF_broadened(TXB_els_HW))*TimeToH;

if nargout>2
    syn_wf=struct('t', t1,'P', TX_WF_broadened, 'mask', TXB_els_HW);
end


% N.B.  Could do this with a LUT as a function of SNR, pulse width, and
% window size.  Would have to be regenerated for each new WF estimate
if 0
    load WF_est;
    sigma_0=0.5*diff(wf_percentile(WF.t, WF.p, [0.16 0.84]));
    W_pulse=sigma_0*1.5e8+[0 0.025 0.05:.05:2];
    
    HW_vals=3+([0:0.1:1, 2:18]);
 
    SNR_vals=logspace(0, 3, 16);
    
    % the minimum SNR is 10 PE /( 10 MHz * 20 m /1.5e8 m/s * 57 pulses) =
    % 1.3
    
    % need a quick test of robust_peak_width_from_hist
    
    %SNR_vals=500;
    for kWp=1:length(W_pulse)
        for kHW=1:length(HW_vals)
            for kSNR=1:length(SNR_vals)
                 [dM(kWp,  kHW, kSNR), dCtr(kWp, kHW, kSNR)]=correct_for_TX_shape(WF.t, WF.p, HW_vals(kHW)/1.5e8,  W_pulse(kWp)/1.5e8, SNR_vals(kSNR));
            end
        end
    end
    TX_shape_corr_table.W_pulse=W_pulse; 
    TX_shape_corr_table.SNR=SNR_vals; 
    TX_shape_corr_table.HW=HW_vals;
    TX_shape_corr_table.delta_med=dM; 
    TX_shape_corr_table.delta_mean=dCtr
end

