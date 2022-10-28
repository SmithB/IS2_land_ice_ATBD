files=glob('/Volumes/ice2/ben/GreenlandSimGF/*v4*/*sigparm*.h5');
fid=fopen('ATL06_queue.txt','w');
for k=1:length(files);
   disp(files{k})
   [out_file, filename]=proc_one_GL_sim_file(files{k}, true);
   if ~exist(out_file,'file');
        fprintf(fid,'proc_one_GL_sim_file(''%s'');\n', files{k});
   end
end	    
