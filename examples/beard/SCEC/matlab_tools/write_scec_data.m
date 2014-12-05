function write_scec_data(filename,header,fields,varargin)

fid = fopen(filename,'w');

for k = 1:length(header)
  fprintf(fid,[header{k},'\n']);
end

fclose(fid);
