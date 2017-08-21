function S=IS_LI_paths(run_type)


% set your data directories here.  Mine are:
S.data_root='/Volumes/ice1/ben/sdt/ATLxx_example';
S.alg_version='/v13_hdf';
S.ATL03=[S.data_root,S.alg_version,'/PIG_ATL03'];
S.ATL06=[S.data_root,S.alg_version,'/PIG_ATL06'];
if exist('run_type','var');
    S.ATL11=[S.data_root,S.alg_version,'/PIG_ATL11_6segs_', run_type];
    S.ATL14=[S.data_root,S.alg_version,'/PIG_ATL14_6segs_', run_type];
end
S.DEM_dir=[S.data_root,'/PIG_DEM'];



