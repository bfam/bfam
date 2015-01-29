function [pdata,data] = process_body_station(base_name, directory)

if nargin < 2
  directory = '.';
end

A = dir(directory);

data(2).name = []; % kills an mlint warning

N_data = 0;
for k = 1:length(A)
  if strncmpi(base_name,A(k).name,length(base_name))
    N_data = N_data + 1;
    data(N_data).name = A(k).name;

    full_name = [directory,'/',A(k).name];

    % Find header lines
    iid = fopen(full_name, 'r');
    header = 0;
    ln = fgetl(iid);
    while ln(1) == '#'
      header = header + 1;
      data(N_data).header{header} = ln;
      ln = fgetl(iid);
    end

    data(N_data).variables = ln;
    if 1 ~= exist('variables','var')
      var_string = ln;
      variables = strsplit(ln,' ');
      num_col = length(variables);
    end
    if strcmp(data(N_data).variables,var_string) ~= 1
      error(['Format different for this file: ',full_name])
    end
    dat = fscanf(iid,'%e');
    data(N_data).dat = reshape(dat,num_col,length(dat)/num_col)';

    fclose(iid);

    for n = 1:length(variables)
      eval(['data(N_data).',variables{n},' = data(N_data).dat(:,n);']);
    end
  end
end
data = data(1:N_data);

if N_data > 0
  n_steps = inf;
  for n = 1:N_data
    n_steps = min(n_steps,length(data(n).t));
  end

  pdata.header = data(1).header;
  for k = 1:length(variables)
    var = ['pdata.',variables{k}];
    eval([var,' = 0;']);
    for n = 1:N_data
      eval([var,' = ',var,'+data(n).dat(:,k)/N_data;']);
    end
  end
else
  pdata = [];
end
