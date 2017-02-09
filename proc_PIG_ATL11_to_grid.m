
% name these conditions
run_type='dh_clouds';

% setup the path defaults
IS_paths=IS_LI_paths('dh_clouds');

 load('/PIG_LaserTracks')

[~, out]=unix(['ls -d ', IS_paths.ATL11, '/run*']);

N_runs=length(strsplit(deblank(out)));


for k_run=1:N_runs
    ATL11_dir=sprintf('%s/run_%d',IS_paths.ATL11, k_run);
    out_file=[ATL11_dir,'/ATL13_solution.mat'];
    
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
    
    
    [M.dx, M.dy]=deal(250);
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
    [grids, G, Gc, TOC, D]=est_dh_dt_fd('build_G',D, M);
    Cvals=est_dh_dt_ft('calc_Cvals', D.sigma, M.sigma, TOC);
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



if false
    [~, out]=unix('ls v9_hdf/PIG_ATL11_dh_clouds/run*/ATL13_solution.mat');
   % [~, out]=unix('ls v9_hdf/PIG_ATL11_with_dh_run_*/tau=2.0//ATL13_solution.mat'); 
    out=strsplit(deblank(out));
    for k=1:length(out); 
        temp=load(out{k},'dzbar'); 
        temp.dzbar.sigma=sqrt(mean((temp.dzbar.z_est(:,:,1:12)-temp.dzbar.z_true(:,:,1:12)).^2,3));
        dzbar(k)=temp.dzbar; 
    end;
    sig=cat(3, dzbar.sigma);
    
    sigma_all=(sqrt(mean(sig.^2,3)));
    
    filename='Dh_real_and_est.gif';
    delete(filename);
    clf; set(gcf,'position', [550 35 1090 1000],'color', [1 1 1],'defaultaxesfontsize', 14,'defaultaxesfontname','times');
    h=cheek_by_jowl(3,1, [0.05 0.07 0.8 0.9]);
    dzg=PIG_dh_function(xg, yg, tg*365);
    dzg=dzg-repmat(dzg(:,:,1), [1 1 size(dzg,3)]);
    dz.z=dz.z-repmat(dz.z(:,:,1), [1 1 size(dz.z,3)]);
    
    MM=double(dz.mask); MM(dz.mask==0)=NaN;
    hb=axes('position', [0.86 0.39 0.02 0.58]); 
    axes(hb); imagesc([0 1], [-11.5:.1:11.5], [-11.5:.1:11.5]'); caxis([-11.5 11.5]);set(hb,'yaxislocation','right','xtick', []); ylabel(hb,'dh WRT h(0), m'); 
    set(gca,'ylim', [-10.5 10.5]);
    colormap([0 0 0; flipud(dzdt_cpt(128))]); 
    set(h(1:2),'visible','off');
    rc1=[14, 42];
    rc2=[14, 27];
    for k=1:13; 
        axes(h(1)); cla
        imagesc(dz.x, dz.y, (dz.z(:,:,k)).*MM); caxis([-11.5 11.5]); axis xy equal tight;  
        XL=get(gca,'xlim'); YL=get(gca,'ylim');
        ht=text(XL*[0.95 0.05]', YL*[0.1 0.9]', sprintf('T=%d months', round(dz.t(k)*365/30.4)),'color', [1 1 1],'fontweight','bold','fontsize', 14);
        hold on;plot(dz.x(rc1(2)), dz.y(rc1(1)),'wo','markersize', 8,'linewidth', 2);
        plot(dz.x(rc2(2)), dz.y(rc2(1)),'wo','markersize', 8, 'linewidth', 2);
        els=abs(D.t-dz.t(k)*365)<.25*365/2;
        
        plot(unique(D.x(els)+1i*D.y(els)),'k.','markersize', 1)
        axes(h(2)); cla;
        imagesc(dz.x, dz.y, (dzg(:,:,k)-dzg(:,:,1)).*MM); caxis([-11.5 11.5]); axis xy equal tight;
        hold on;plot(dz.x(rc1(2)), dz.y(rc1(1)),'wo','markersize', 8,'linewidth', 2);
        plot(dz.x(rc2(2)), dz.y(rc2(1)),'wo','markersize', 8, 'linewidth', 2);
        plot(unique(D.x(els)+1i*D.y(els)),'k.','markersize', 1)
        set(h(1:2),'visible','off');
        set(findobj(gca,'type','line','marker','o'),'markerfacecolor', [0 0 0]);
        
        
        axes(h(3)); cla; hold on
        plot( (dz.t*365/30.4), squeeze(dzg(rc1(1), rc1(2), :)),'ko'); 
        plot( (dz.t(1:k)*365/30.4), squeeze(dz.z(rc1(1), rc1(2), 1:k)),'b','linewidth', 3);
        plot( (dz.t*365/30.4), squeeze(dzg(rc2(1), rc2(2), :)),'ko'); 
        plot( (dz.t(1:k)*365/30.4), squeeze(dz.z(rc2(1), rc2(2), 1:k)),'r','linewidth', 3);
        ylabel('dh WRT h(0), m');xlabel('month');
        set(gca,'xlim', [-0.5 36.5] ,'ylim', [-10.5 10.5]); 
        grid on;
        
        drawnow; frame=getframe(gcf); 
        im=frame2im(frame);
        [A, map]=rgb2ind(im, 256); 
        if k==1; 
            imwrite(A, map, filename,'gif', 'LoopCount', Inf,'DelayTime', 1);
        else
            imwrite(A, map, filename,'gif', 'DelayTime', 1,'WriteMode', 'append');
        end
    end
  
    imwrite(A, map, filename,'gif', 'DelayTime', 5,'WriteMode', 'append');
    
    
    
    % spin the DEMs movie;
    II=read_geotif('PIG_DEM/20140119_1453_10200100282DE200_102001002A716000-DEM_tr16x.tif');
    II.z(II.z==0)=NaN;
    M0=double(z0.mask); M0(z0.mask==0)=NaN;
    ha=axes('position', [0 0 1 1]);
    hold on;
    hs(1)=surf(II.x(1:10:end), II.y(1:10:end), double(II.z(1:10:end, 1:10:end)));shading interp; material dull;  hs(2)==surf(z0.x, z0.y-3e4, full(z0.z.*M0)); shading interp; material dull; lighting phong;  hl=light;
    set(gca,'dataaspectratio', [ 1 1 .01]);
    set(gca,'visible','off')
    view(0, 70);
    filename='two_DEM_movie.gif';
    for theta=0:2:360;
        view( theta, 70);
        drawnow; frame=getframe(gcf); 
        im=frame2im(frame);
        [A, map]=rgb2ind(im, 256); 
        if theta==0; 
            imwrite(A, map, filename,'gif', 'LoopCount', Inf,'DelayTime', 1/8);
        else
            imwrite(A, map, filename,'gif', 'DelayTime', 1/8,'WriteMode', 'append');
        end
    end
  
    imwrite(A, map, filename,'gif', 'DelayTime', 5,'WriteMode', 'append');
    
    
    % make an error map
    [II.gx]=gradient(II.z, II.x, II.y);
    figure(123);clf
    imagesc(II.x(1:5:end), II.y(1:5:end), repmat(scale_to_byte(II.gx(1:5:end, 1:5:end), [-.05 0.05]), [1 1 3]),'alphadata', double(isfinite(II.gx(1:5:end, 1:5:end)))); 
    axis xy equal tight
    hold on;
    plot(unique(D.x+1i*D.y),'k.','markersize', 1);
    imagesc(dzbar(1).x, dzbar(1).y, sigma_all,'alphadata', 0.7*isfinite(sigma_all));
    hb=colorbar; 
    ylabel(hb,'RMSE(\deltaz), m');
    set(gca,'visible','off');
    figure(124)
    hold on;
    for k=1:10; plot(dz.t, squeeze(dzbar(k).z_est(2, 3,:))); end
    plot(dz.t, squeeze(dzbar(1).z_true(2,3,:)),'ko'); 
    xlabel('year'); ylabel('\delta h, m'); 
    
    
end