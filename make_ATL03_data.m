%-------------------------------------------------------
function [D2, params]=make_ATL03_data(N_pulses, N_chan, roughness, N_per_pulse, WF, BGR, H_window, AT_slope)

if ~exist('H_window','var');
    H_window=4;
end

if exist('AT_slope','var');
    DEM=struct('x', -50:50:(50*ceil(N_pulses*.7/50+1)), 'y', [-100 100]');
    [xg,~]=meshgrid(DEM.x, DEM.y);
    DEM.z=xg*AT_slope;
else
    DEM=0;
end
x0=(1:N_pulses)*0.7;
ATM_xmit=ones(size(x0));
params=struct('roughness', roughness,   'sigma_x', 7.5, 'NoiseRate', BGR, 'H_window', H_window, 'WF', WF, 'c', 3e8, 'N_channels', N_chan, 'N_per_pulse', N_per_pulse,'t_dead', 3.2e-9);
D2=det_sim(DEM, x0, params);
