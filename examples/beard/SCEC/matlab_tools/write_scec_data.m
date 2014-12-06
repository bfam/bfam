function write_scec_data(filename,header,fields,DAT)

fid = fopen(filename,'w');

for k = 1:length(header)
  fprintf(fid,[header{k},'\n']);
end
fprintf(fid,[fields,'\n']);

for k = 1:size(DAT,1)
  for n = 1:size(DAT,2)
    fprintf(fid,'%+30.17E',DAT(k,n));
    if(n ~= size(DAT,2))
      fprintf(fid,' ');
    end
  end
  fprintf(fid,'\n');
end

fclose(fid);
