
% name these conditions
run_type='dh_clouds';

% setup the path defaults
IS_paths=IS_LI_paths('dh_clouds');

 load('/PIG_LaserTracks')

[~, out]=unix(['ls -d ', IS_paths.ATL11, '/run*']);

N_runs=length(strsplit(deblank(out)));


for k_run=1:N_runs
    ATL11_dir=sprintf('%s/run_%d',IS_paths.ATL11, k_run);

    ATL14_dir=sprintf('%s/run_%d',IS_paths.ATL14, k_run);
    if ~exist(ATL14_dir,'dir'); mkdir(ATL14_dir); end
    out_file=[ATL14_dir,'/ATL14_solution.mat'];
    
    % !!! skip the lock check !!!
    %if lockfile_tool('lock', out_file); continue; end;
    %if exist(out_file,'file'); continue; end;
    
    % loop over data files for the current run
    fprintf(1, 'working on %s\n', out_file);
    [s, files]=unix(sprintf('ls %s/*D3b.h5', ATL11_dir));
    files=strsplit(deblank(files));
    clear D11 tau track pair
    [track, pair]=deal(NaN(length(files), 1));
    % load the ATL11 data
    for k=1:length(files)
        % parse the segment and track from the file name

        temp=regexp(files{k}, 'Track_(\d+)-Pair_(\d+)_\S+.h5','tokens');
        track(k)=str2num(temp{1}{1});
        pair(k)=str2num(temp{1}{2});
        
        temp=read_ATL11_h5(files{k}); 
        
        % convert the lat and lon measurement positions to x and y
        xx=ll2ps(temp.corrected_h.ref_pt_lat, temp.corrected_h.ref_pt_lon);
        temp.corrected_h.ref_pt_x_PS=real(xx);
        temp.corrected_h.ref_pt_y_PS=imag(xx);
         
        h1=temp.corrected_h.pass_h_shapecorr;
        s1=temp.corrected_h.pass_h_shapecorr_sigma;
        mask=isfinite(h1) & isfinite(s1);
        W=1./s1.^2;
        W(s1==0 | mask==0)=0;
        hbar=sum(h1.*W, 2)./sum(W, 2);
        hbar(sum(W,2)==0)=NaN;
        temp.corrected_h.dz=h1-repmat(hbar, [1, 12]);
        temp.reference_point.PT=ones(size(temp.reference_point.seg_count))*pair(k);
        temp.reference_point.RGT=ones(size(temp.reference_point.seg_count))*track(k);
        D11(k)=temp;
    end
    
    clear Dcorr
    temp=cat(1, D11.corrected_h);
    ff=fieldnames(temp);
    for kf=1:length(ff)
        Dcorr.(ff{kf})=cat(1, temp.(ff{kf}));
    end
    temp=cat(1, D11.reference_point);
    Dcorr.PT=cat(1, temp.PT);
    Dcorr.RGT=cat(1, temp.RGT);
    
    
    clear D;
    D.x=repmat(Dcorr.ref_pt_x_PS,  [1, 12]);
    D.y=repmat(Dcorr.ref_pt_y_PS,  [1, 12]);
    D.rep=repmat(1:12, [size(Dcorr.ref_pt_x_PS,1),1]);
    D.x=real(D.x);
    D.t=Dcorr.mean_pass_time;
    D.h=Dcorr.pass_h_shapecorr;
    D.sigma=sqrt(Dcorr.pass_h_shapecorr_sigma.^2+Dcorr.pass_h_shapecorr_sigma_systematic.^2);
    D.PT=repmat(Dcorr.PT, [1 12]);
    D.RGT=repmat(Dcorr.RGT, [1 12]);
  
    
    D=index_struct(D, isfinite(D.x) & isfinite(D.h) & isfinite(D.sigma) & abs(D.sigma) < 2);
    
    
    ff=fieldnames(D);
    for kf=1:length(ff)
        D.(ff{kf})=D.(ff{kf})(:);
    end
    
    
    [M.dx, M.dy]=deal(500);
    [M.dx0, M.dy0]=deal(125);
    M.dt=0.25;
    
    M.XR=round_to(range(D.x), M.dx)+[-M.dx M.dx];
    M.YR=round_to(range(D.y), M.dx)+[-M.dx M.dx];
    M.TR=[0 3.25];
    M.time_zero_season=6;
    
    W=sqrt(diff(M.YR).*diff(M.XR));
    
    params=struct('W', W,'dx', M.dx, 'dx0', M.dx0, 'dt', M.dt, ...
        'est_errors', false,'time_zero_season', 6);%, 'calc_PS_bias', true);
    %use default params for selected for Thwaties lakes: probably lax for the d2zdt2 parameter.  Set these as the
    %default
    params.E_RMS_d2z0dx2=2e-7;
    params.E_RMS_d3zdx2dt=6e-8;
    params.E_RMS_d2zdt2=1;
    
    % define sigma structure:
    line_spacing_z0=1e3;
    line_spacing_dz=1e3;
    sigma=struct( 'smooth_dz',params.E_RMS_d3zdx2dt*W, 'smooth_z0', params.E_RMS_d2z0dx2*W, 'same_season_dz', 2*W,...
        'smooth_season_dz', params.E_RMS_d2zdt2*W,'smooth_bias', .001*W, 'zero_bias', .001*W, 'line_spacing_z0', line_spacing_z0, 'line_spacing_dz', line_spacing_dz );
    sigma.time_gap=0.25;
        
    M.sigma=sigma;
    
    D.year=D.t/365;
    [G, Gc, TOC, grids, D]=est_dh_dt_fd('build_G',D, M);
    Cvals=est_dh_dt_fd('calc_Cvals', D.sigma, M.sigma, TOC);
    d_c=zeros(size(Gc(:,1)));
    good=true(size(D.h));
    
    % iterate to remove outliers.  five iterations is arbitrary
    for k=1:5
        C=spdiags([Cvals.G(good,:);  Cvals.Gc], 0, sum(good)+size(Gc,1), sum(good)+size(Gc,1));
        % run the fit
        tic;
        m=my_lscov([G(good,:); Gc], double([D.h(good); d_c]), C);this_dt=toc;
        r=D.h-G*m;
        this_sigma=iqr(r./D.sigma)/2;
        % eliminate data that are more than three sigma from the model
        good=abs(r./D.sigma)< 3*this_sigma & abs(r) < 5;
        fprintf(1, 'k=%d, sigma=%3.2f; culling %3.2f percent\n', k, this_sigma, 100*(1-mean(good)));
    end
    
    % calculate the residuals
    D.d_est=G*m;
    D.r=D.h-D.d_est;
    D.good=good;
    
    
    % build the output grid
    dz.y=grids.dz.ctrs{1};
    dz.x=grids.dz.ctrs{2};
    dz.t=grids.dz.ctrs{3};
    dz.z=reshape(m(TOC.cols.dz), grids.dz.dims);
    
    % mask of constrained dz nodes (the rest are filled by interpolated
    % values)
    [xg, yg]=meshgrid(dz.x, dz.y);
    mask=false(size(xg));
    mask(ismember(xg+1i*yg, round_to(D.x+1i*D.y, M.dx)))=true;
    dz.mask=dilate_mask(mask, 4);
    
    % mask of constrained DEM points
    z0.x=grids.z0.ctrs{2};
    z0.y=grids.z0.ctrs{1};
    z0.z=reshape(m(TOC.cols.z0), grids.z0.dims);
    [xg0, yg0]=meshgrid(z0.x, z0.y);
    mask=false(size(xg0));
    mask(ismember(xg0+1i*yg0, round_to(D.x+1i*D.y, M.dx)))=true;
    z0.mask=double(dilate_mask(mask, 4*floor(M.dx/M.dx0)));
    z0.mask(z0.mask==0)=NaN;
    
    % display the output
    hs=surf(z0.x, z0.y, z0.z.*z0.mask); shading interp; set(gca,'dataaspectratio', [1 1 .01]);material dull;  hl=light;
    
    % evaluate the output: Generate the expected elevation-change function
    [yg,xg,tg]=ndgrid(dz.y, dz.x, dz.t);
    dzg=PIG_dh_function(xg, yg, tg*365);
    dzg=dzg-repmat(dzg(:,:,6), [1 1 size(dzg,3)]);
    
    % we only need one layer of xg and yg
    xg=xg(:,:,1); yg=yg(:,:,1);
    
    % build a structure called dzbar
    clear dzbar
    [xg1, yg1]=meshgrid(min(xg(1,:)):5000:max(xg(1,:)), min(yg(:,1)):5000:max(yg(:,1)));
    dzbar.x=xg1(1,:);
    dzbar.y=yg1(:,1);
    dzbar.count=zeros(size(xg1))+NaN;
    dzbar.z_true=NaN([size(xg1),size(tg,3)]);
    dzbar.z_est=dzbar.z_true;
    for k=1:length(xg1(:))
        els=abs(xg+1i*yg- (xg1(k)+1i*yg1(k)))<=5700;
        [r,c]=ind2sub(size(xg1(:,:,1)), k);
        these=dz.mask & els;
        dzbar.count(r,c)=sum(these(:));
        if sum(these(:))>90
            for kt=1:size(tg,3)
                temp=dzg(:,:,kt);
                temp1=dz.z(:,:,kt);
                dzbar.z_est(r, c, kt)=mean(mean(temp1(these)));
                dzbar.z_true(r, c, kt)=mean(mean(temp(these)));
            end
        end
    end
    
    
    dzbar.sigma=sqrt(mean((dzbar.z_est-dzbar.z_true).^2,3));
    
    save( out_file, 'r','good','D','m','TOC','M', 'grids','dz','z0','dzbar');
    lockfile_tool('unlock', out_file); 
end


 