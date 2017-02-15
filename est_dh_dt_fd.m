function varargout=est_dh_dt_fd(varargin)

% varargout=est_dh_dt_fd(varargin)
% wrapper containing functions to estimate a time-zero DEM (z0) and a set
% of elevation-change surfaces (dz).  These functions get called by 
% proc_PIG_ATL11_to_grid and est_cryosat_dzdt


if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end


%-----------------------------------------------------------------
function  [G, Gc,  TOC, grids, D]=build_G(D, M)

% [G, Gc,  TOC, grids, D]=build_G(D, M)
% build the least-squares fitting matrix based on the data and the grid-
% definition structure M.
% inputs:
% D:  input data structure.  Must have, at minimum, fields:
%        x, y, h, t, sigma
% M:  grid-definition structure, containing, at minimum fields:
%    dx, dy: horizontal spacing for elevation change nodes
%    dx0, dy0: "          "      "  DEM nodes
%    dt:    temporal spacing of elevation change surfaces
%    XR:    range of x values
%    YR:    range of y values
%    TR:    time range 
%    time_zero_season:  time of reference DEM
 
% outputs:
% grids:  Structure defining the output grids
% G:      Matrix that, multiplied by the model vecotor, provides
%         reconstructed data values
% Gc:     Structure giving constraint matrices that, multipled by the model vector,  gives the derivatives of the model
% TOC:    Structure defining the equation types in Gc
% D:      Input data, subset to contain points that fall within the model
%         domain

% note that the constraint equations here are normalized so that for each
% component, if you add the squares of G*m over all nodes, you get the spatial
% integral of the square of the quantity.

eq_key=define_equation_types;

%dz grid
nx=ceil(diff(range(M.XR))./M.dx)+1;
ny=ceil(diff(range(M.YR))./M.dy)+1;
nt=ceil(diff(range(M.TR))./M.dt)+1;
grids.dz.dims=[ny nx nt];
Ndz=nx*ny*nt;

grids.dz.ctrs{1}=M.YR(1)+(0:ny-1)*M.dy;
grids.dz.ctrs{2}=M.XR(1)+(0:nx-1)*M.dx;
grids.dz.ctrs{3}=M.TR(1)+(0:nt-1)*M.dt;
TOC_cols.dz=1:prod(grids.dz.dims);
TOC_cols.dz_epoch=repmat(1:nt, [ny*nx, 1]);
TOC_cols.dz_epoch=TOC_cols.dz_epoch(:)';

%z0 grid
nx0=ceil(diff(range(M.XR))./M.dx0)+1;
ny0=ceil(diff(range(M.YR))./M.dy0)+1;
grids.z0.dims=[ny0 nx0];
grids.z0.ctrs{1}=M.YR(1)+(0:ny0-1)*M.dy0;
grids.z0.ctrs{2}=M.XR(1)+(0:nx0-1)*M.dx0;
TOC_cols.z0=Ndz+(1:prod(grids.z0.dims));

Ncols=prod(grids.dz.dims)+prod(grids.z0.dims);

% delete D values for which the interpolations will be undefined.
good=true(size(D.x));
good=good & D.x>=grids.dz.ctrs{2}(1) & D.x <= grids.dz.ctrs{2}(end);
good=good & D.y>=grids.dz.ctrs{1}(1) & D.y <= grids.dz.ctrs{1}(end);
good=good & D.year>=grids.dz.ctrs{3}(1) & D.year <= grids.dz.ctrs{3}(end);
D=index_struct(D, good);

k=0;
%smooth dz/dt maps
k=k+1;
dx=M.dx;
dt=M.dt;
d3zdx2dt.sub={[0 0 0 0 0 0], [-1 0 1 -1 0 1], [-1 -1 -1 0 0 0]};
d3zdx2dt.val=[-1 2 -1 1 -2 1]/M.dx/dt; % integrand
CEQ(k)=make_fd_op(grids.dz, d3zdx2dt, 0,  021, eq_key.smooth_dz);
k=k+1;
d3zdy2dt.sub={[-1 0 1 -1 0 1], [0 0 0 0 0 0], [-1 -1 -1 0 0 0]};
d3zdy2dt.val=[-1 2 -1 1 -2 1]/M.dx/dt; % integrand
CEQ(k)=make_fd_op(grids.dz,  d3zdy2dt, 0, 201, eq_key.smooth_dz);
k=k+1;
d3zdxdydt.sub={[0, 1, 0, 1, 0, 1, 0, 1], [0, 0, 1, 1, 0, 0, 1, 1], [-1 -1 -1 -1 0 0 0 0]};
d3zdxdydt.val=[-1 1 1 -1 1 -1 -1 1]/M.dx/M.dt; % integrand
CEQ(k)=make_fd_op(grids.dz, d3zdxdydt, 0, 111, eq_key.smooth_dz);

% flat dz/dt maps
k=k+1;
d2zdxdt.sub={ [0 0 0 0], [-1 0 -1 0], [-1 -1 0 0]};
d2zdxdt.val=[-1 1 1 -1]/dt; % integrand
CEQ(k)=make_fd_op(grids.dz, d2zdxdt, 0, 11, eq_key.flat_dz);
k=k+1;
d2zdydt.sub={[-1 0 -1 0], [0 0 0 0],  [-1 -1 0 0]};
d2zdydt.val=[-1 1 1 -1]/dt; % integrand
CEQ(k)=make_fd_op(grids.dz, d2zdydt, 0,  101,  eq_key.flat_dz);

% smooth dz in time
k=k+1;
d2zdxdt.sub={ [0 0 0], [ 0 0 0], [-1 0 1]};
d2zdxdt.val=[-1 2 -1]/dt/dt*dx; % integrand
CEQ(k)=make_fd_op(grids.dz, d2zdxdt, 0, 2, eq_key.smooth_season_dz);

if isfield(M.sigma,'time_gap');
    % flat dz in time
    k=k+1;
    d2zdxdt.sub={ [0 0 ], [ 0 0 ], [-1 1]};
    d2zdxdt.val=[-1 1]/dt*dx; % integrand
    CEQ(k)=make_fd_op(grids.dz, d2zdxdt, 0, 1, eq_key.flat_season_dz);
end

% zero dz at time zero
k=k+1;
CEQ(k).r=(1:nx*ny)';
CEQ(k).c=find(TOC_cols.dz_epoch==M.time_zero_season)';
CEQ(k).v=ones(nx*ny,1);
CEQ(k).eqn_id=eq_key.zero_dz_at_t0*ones(nx*ny,1);
CEQ(k).eqn_type=eq_key.zero_dz_at_t0*ones(nx*ny,1);
CEQ(k).node_num=(1:nx*ny)';
 
% smooth z0
k=1;
dx0=M.dx0;
d2z0dx2.sub={[0 0 0 ], [-1 0 1]};
d2z0dx2.val=[-1 2 -1]/dx0;   % N.B. This is the sqare root of the spatial integrand. This would be /dx0/dx0, but the spatial integral cancels out one of the deltas
CEQ0(k)=make_fd_op(grids.z0, d2z0dx2,  Ndz, 020, eq_key.smooth_z0);
k=k+1;
d2z0dy2.sub={ [-1 0 1], [0 0 0]};
d2z0dy2.val=[-1 2 -1]/dx0; % integrand
CEQ0(k)=make_fd_op(grids.z0, d2z0dy2,  Ndz, 200, eq_key.smooth_z0);

k=k+1;
d2z0dxdy.sub={[0 1 0 1], [0 0 1 1 ]};
d2z0dxdy.val=[-1 1 1 -1]/dx0; % integrand
CEQ0(k)=make_fd_op(grids.z0, d2z0dxdy,  Ndz, 110, eq_key.smooth_z0);

% flat z0
k=k+1;
dz0dx.sub={[0 0], [-1 1]};
dz0dx.val=[-1 1]; % integrand
CEQ0(k)=make_fd_op(grids.z0, dz0dx,  Ndz, 010, eq_key.flat_z0);
k=k+1;
dz0dy.sub={ [-1 1], [0 0]};
dz0dy.val=[-1 1]; % integrand
CEQ0(k)=make_fd_op(grids.z0, dz0dy,  Ndz, 100, eq_key.flat_z0);

%Fitting matrices for dz
A_interp_dz=bilinear_basis_functions(double(D.x), double(D.y),  grids.dz.ctrs{2}, grids.dz.ctrs{1});
W=bilinear_basis_functions(double(D.year), zeros(size(D.year)), grids.dz.ctrs{3}(:), 0);
%W=fd_spline_fit('build_G',  D.year, grids.dz.ctrs{3}, 0);

% build the fitting matrix for dz
% structure of G is  [t1:dz_1...dz_n], [t2:dz_1...dz_n], [t3:dz_1...dz_n], etc
% then the column for each point is c_i + (c_t-1)*Nnodes
for k=1:nt
    % multiply each column of A_interp by the coefficients of the time interpolation function for year k

    if datenum(version('-date')) > datenum('June 1 2016')
        % newer versions of Matlab use the broadcasting rule to allow
        % multiplication of a vector with an array.  This is fast and
        % efficient in this case
        G0=A_interp_dz.*W(:,k);
    else
        % OTW have to generate a replicated version of W, which is slow
        G0=A_interp_dz.*repmat(W(:,k), [1, size(A_interp_dz, 2)]);
    end
    [Gdz(k).r,Gdz(k).c,Gdz(k).v]=find(G0);
    Gdz(k).c=Gdz(k).c+nx*ny*(k-1);
    %toc
end
Gdz=sparse(cat(1, Gdz(:).r), cat(1, Gdz(:).c), cat(1, Gdz(:).v), length(D.x), Ncols);
 
% fitting matrix for z0
[Gz0.r, Gz0.c, Gz0.v]=find( bilinear_basis_functions(D.x, D.y,  grids.z0.ctrs{2}, grids.z0.ctrs{1}));
Gz0.c=Gz0.c+prod(grids.dz.dims);
Gz0.eqn_id=3+zeros(length(D.x),1);
Gz0.eqn_type=eq_key.fit+zeros(length(D.x(:)),1);

% assemble the z0 and dz fitting matrices
G=Gdz+sparse(Gz0.r, Gz0.c, Gz0.v, length(D.x(:)), Ncols);
%[TOC.eqn_id, TOC.eqn_type]=deal(zeros(size(G,1),1)+eq_key.fit);
  
[Gc, TOC]=assemble_vert([CEQ, CEQ0], Ndz+prod(grids.z0.dims));
TOC.cols=TOC_cols;




%-----------------------------------------------------------------
function Cvals=calc_Cvals(data_sigma, sigma, TOC)

% Cvals=prep_fit(data_sigma, sigma)
% Return the diagonal to the covariance matrix for the model fit
% inputs:
%   data_sigma: Ndatax1 vector of error estimates for the data values
%   sigma: sturcture containing constraint-equaiton scales, with fields:
%         smooth_dz:   expected square integral per node of d2(dz)/dtd(x or y)2
%         line_spacing_dz:   Size of data gaps in dz (line to line)
%         smooth_z0: expected square integral per node of d2(dz)/d(x or y)2
%         line_spacing_z0: Size of data gaps in z0 (line to line)
%         smooth_season_dz: expected square integral per node of d2(dz)/dt2
%         and, optionally: 
%             time_gap: the expected time gap between data 
% outputs:
% Cvals:  Structure giving expected values for the squared output values


eq_key=define_equation_types;
Cvals.G=data_sigma(:).^2;
Cvals.Gc=zeros(size(TOC.eqn_type));
Cvals.Gc(TOC.eqn_type==eq_key.smooth_dz)=sigma.smooth_dz.^2;
Cvals.Gc(TOC.eqn_type==eq_key.flat_dz)=sigma.smooth_dz.^2*sigma.line_spacing_dz.^2;
Cvals.Gc(TOC.eqn_type==eq_key.smooth_z0)=sigma.smooth_z0.^2;
Cvals.Gc(TOC.eqn_type==eq_key.flat_z0)=sigma.smooth_z0.^2*sigma.line_spacing_z0.^2;
Cvals.Gc(TOC.eqn_type==eq_key.zero_dz_at_t0)=sigma.smooth_z0.^2/1e4;
Cvals.Gc(TOC.eqn_type==eq_key.smooth_season_dz)=sigma.smooth_season_dz.^2;
if isfield(sigma,'time_gap')
    Cvals.Gc(TOC.eqn_type==eq_key.flat_season_dz)=sigma.smooth_season_dz.^2*sigma.time_gap.^2;
end
 

%-----------------------------------------------------------------
function Gc_plots(S, which_plots)

% make plots showing the magnitude of different components of the
% constraint equations for a model structure

if ~exist('which_plots','var')
    which_plots=20;
end

grids=S.grids;
figure; clf;
surfl(grids.z0.ctrs{1}, grids.z0.ctrs{1}, S.z0);
N_seas=size(S.dz,3);
N=ceil(sqrt(size(S.dz,3)));
figure; clf
h=cheek_by_jowl(N,N,[0 0 0.9 0.9]);
for k=1:N_seas; axes(h(k)); imagesc(S.dz(:,:,k)); caxis([-1 1]); end

if which_plots<3;
    return
end
f=fieldnames(S.ceq);
for kf=1:length(f);  
     for kk=1:length(S.ceq.(f{kf}));
         figure
        sz=size(S.ceq.(f{kf})(kk).r);
        if length(sz)==2
            n=1;
            sz(end+1)=1;
        else
            n=ceil(sqrt(sz(3)));
        end
        h=cheek_by_jowl(n, n, [0 0 0.9 0.9]);
        for k=1:sz(3);
            axes(h(k));
            imagesc(S.ceq.(f{kf})(kk).r(:,:,k));
            caxis(range(S.ceq.(f{kf})(kk).r(:))); 
        end
        if k < sz(3); 
            delete(h(k+1:end)); 
        end
        colorbar
        axes(h(1)); title([f{kf}, ' ',num2str(S.ceq.(f{kf})(kk).eq_id)]);
    end
end

 
%-----------------------------------------------------------------
function eq_key=define_equation_types

% The constraint equations come with a TOC (table of contents) vector
% defining what kind of equation is found in each row.  

eq_key=struct('fit', 1, ...
    'smooth_dz', 2, ...
    'flat_dz', 2.1, ...
    'smooth_z0', 3, ...
    'flat_z0', 3.1, ...
    'zero_dz_at_t0', 4, ...
     'smooth_season_dz', 6, ...
     'flat_season_dz', 6.1, ...
     'smooth_bias', 7, ...
     'flat_bias', 7.1, ...
     'zero_bias', 7.2, ...
     'seasonal_cycle', 8);
    %'same_season_dz', 5, ...
  

%-----------------------------------------------------------------
function R2=parse_residuals(r, Cvals, TOC);

% Take the RSSS of each component of the residual vector-- the extra S 
% is 'scaled'

eq_key=define_equation_types;
f=fieldnames(eq_key);
rs=r./sqrt(Cvals);
for k=1:length(f);
    R2.(f{k})=sum(rs(TOC.eqn_type==eq_key.(f{k})).^2);
end
R2.dz=R2.smooth_dz+R2.flat_dz;
R2.z0=R2.smooth_z0+R2.flat_z0;
R2.total=R2.dz+R2.z0+R2.smooth_season_dz;%+R2.same_season_dz

%-----------------------------------------------------------------
function S=parse_fit(G, G_c, m, d, TOC, grids)

% Provide grids mapping each component of the model.  In the fitting 
% process, the model grids are flattened into a model vector; this function
% reverses the processs.


S.dz=reshape(m(TOC.cols.dz), grids.dz.dims);
S.z0=reshape(m(TOC.cols.z0), grids.z0.dims);

S.d_est=G*m;
S.d=d;
S.r_d=d-S.d_est;
eq_key=define_equation_types;
f=fieldnames(eq_key);

r_c=G_c*m;

for k=1:length(f);
    if ~isempty(findstr(f{k},'dz'));
        shp=grids.dz.dims;
    elseif ~isempty(findstr(f{k},'z0'));
        shp=grids.z0.dims;
    elseif strcmp(f{k},'fit') || strcmp(f{k},'seasonal_cycle') || strcmp(f{k},'roll_fit') || strcmp(f{k},'roll_fit') || strcmp(f{k},'power_fit');
        continue
    end
    these=TOC.eqn_type==eq_key.(f{k});
    eq_ids=unique(TOC.eqn_id(these));
    if ~any(these); continue; end
    eqn_ids=unique(TOC.eqn_id(these));
    for k_id=1:length(eqn_ids)
        these_id=these & TOC.eqn_id==eqn_ids(k_id);
        S.ceq.(f{k})(k_id).r=zeros(shp);
        S.ceq.(f{k})(k_id).r(TOC.node_num(these_id))=r_c(these_id);
        S.ceq.(f{k})(k_id).eq_id=eqn_ids(k_id);
    end
    S.r.(f{k})=r_c(these);
    S.eqn_id.(f{k})=TOC.eqn_id(these);
end
 %-------------------------------------------------------------------------
function [M, TOC]=assemble_vert(RCV, Ncols)
% given arrays in r, c, v form, stick them together vertically
last_row=0;
row_list=cell(length(RCV),1);
for k=1:length(RCV);
    RCV(k).r=RCV(k).r+last_row;
    row_list{k}=(last_row+1):max(RCV(k).r);
    last_row=max(RCV(k).r);
end
M=sparse(cat(1,RCV.r), cat(1,RCV.c), cat(1,RCV.v), last_row, Ncols);
TOC.eqn_id=cat(1, RCV(:).eqn_id);
TOC.eqn_type=cat(1, RCV(:).eqn_type);
if isfield(RCV,'node_num');
    TOC.node_num=cat(1, RCV(:).node_num);
end

%-------------------------------------------------------------------------
function [G_fit_bias, Gc_bias, grids, TOC_bias]=bias_interp_mtx(M, D, weight, eq_key)

%z0 grid
nx=ceil(diff(range(M.XR))./M.dx)+1;
ny=ceil(diff(range(M.YR))./M.dy)+1;
grids.bias.dims=[ny nx];
grids.bias.ctrs{1}=M.YR(1)+(0:ny-1)*M.dy;
grids.bias.ctrs{2}=M.XR(1)+(0:nx-1)*M.dx;

G_fit_bias=bilinear_basis_functions(D.x, D.y,  grids.bias.ctrs{2}, grids.bias.ctrs{1});
for k=1:size(G_fit_bias,2);
    G_fit_bias(:,k)=G_fit_bias(:,k).*weight;
end

k=1;
% smooth bias
dx=M.dx;
d2zdx2.sub={[0 0 0 ], [-1 0 1]};
d2zdx2.val=[-1 2 -1]/dx; % integrand
CEQ(k)=make_fd_op(grids.bias, d2zdx2, 0, 020, eq_key.smooth_bias);
k=k+1;
d2zdy2.sub={ [-1 0 1], [0 0 0]};
d2zdy2.val=[-1 2 -1]/dx; % integrand
CEQ(k)=make_fd_op(grids.bias, d2zdy2, 0, 200, eq_key.smooth_bias);
k=k+1;
d2zdxdy.sub={[0 1 0 1], [0 0 1 1 ]};
d2zdxdy.val=[-1 1 1 -1]/M.dx; % integrand
CEQ(k)=make_fd_op(grids.bias, d2zdxdy, 0, 110, eq_key.smooth_bias);

% flat bias
k=k+1;
dzdx.sub={[0 0], [-1 1]};
dzdx.val=[-1 1]; % integrand
CEQ(k)=make_fd_op(grids.bias, dzdx, 0, 010, eq_key.flat_bias);
k=k+1;
dzdy.sub={ [-1 1], [0 0]};
dzdy.val=[-1 1]; % integrand
CEQ(k)=make_fd_op(grids.bias, dzdy, 0, 100, eq_key.flat_bias);

k=k+1; 
b_zero.sub={0, 0};
b_zero.val=dx;
CEQ(k)=make_fd_op(grids.bias, b_zero, 0, 0, eq_key.zero_bias);

[Gc_bias, TOC_bias]=assemble_vert(CEQ, prod(grids.bias.dims));

%-------------------------------------------------------------------------
function CEQ=make_fd_op(grids, template, col0, eq_id, eqn_type)

% CEQ=make_fd_op(grids, template, col0, eq_id, eqn_type)
% builds a constraint equation based on the discrete representation of a 
% differential operator.  The user provies a template, giving the node
% offsets around the central node and the coefficients applied to the
% values at each offset.  The returned sparse matrix, when multiplied by
% a flattened model vector, gives the operator value at each node.  
% Boundary nodes are not treated- if the template extends over the edges of
% the grid, no equation is returned.

d_ind=cumprod([1, grids.dims]);
ind=(1:prod(grids.dims))';
Nnodes=ind(end);

ind1=repmat(ind, [1, numel(template.val)]);
row=repmat(ind, [1, numel(template.val)]);
for k=1:length(template.sub)
    ind1=ind1+d_ind(k)*repmat(template.sub{k}, [Nnodes,1]);
end

% put the deltas together
delta_sub=cat(1, template.sub);

node_num=(1:Nnodes)';
% look for bad offsets
bad=any(ind1<1 | ind1 > Nnodes,2);
row=row(~bad,:);
ind1=ind1(~bad,:);
node_num=node_num(~bad,:);
% look for wrapped offsets.  Get the offset by taking the difference between 
% row and ind1 and converting it to an index based on the shape of the offset matrix.
% This should equal the offset implied in the template.
ind_0=cell(length(grids.dims), 1);
[ind_0{:}]=ind2sub(grids.dims, row);
ind_off=cell(length(grids.dims), 1);
[ind_off{:}]=ind2sub(grids.dims, ind1);
delta_ind=cat(3, ind_off{:})-cat(3, ind_0{:});
bad=any(any(delta_ind ~= permute(repmat(cat(1, delta_sub{:}), [1 1 size(ind1,1)]), [3 2 1]), 2),3);

ind1=ind1(~bad,:);
N_rows=sum(~bad);
row=row(~bad,:);
node_num=node_num(~bad);
% eliminate blank rows.
[ur, ~, ind]=unique(row);
row=ind;
val=repmat(template.val, [sum(~bad), 1]);
CEQ.c=ind1(:)+col0; 
CEQ.v=val(:);
CEQ.r=row(:);
CEQ.eqn_id=eq_id*ones(N_rows, 1);
CEQ.eqn_type=eqn_type*ones(N_rows, 1);
CEQ.node_num=node_num(:);


 

  