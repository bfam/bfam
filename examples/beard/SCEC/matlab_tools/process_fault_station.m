% n_ref is the positive normal direction. All that is used is the sign of
% n_ref \cdot n and using this, if n = [0,1,0] and n_ref=[0,1,0]:
%      vP1 = vxP -  vxM
%      Vp2 = 0 (no opening)
%      vP3 = vzP -  vzM
function [pdata,data] = process_fault_station(n_ref,base_name, directory)

if nargin < 2
  directory = '.';
end

A = dir(directory);

data(2).name = []; % kills an mlint warning

variables = {'t','Tp1','Tp2','Tp3','Tn','V','Vp1','Vp2','Vp3','Dp','Dp1','Dp2','Dp3','Dn'};

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
      if strcmp(ln,'# interpolated normal:') == 1
        header = header + 1;
        data(N_data).header{header} = ln;
        ln = fgetl(iid);
        data(N_data).n = sscanf(ln,'#    (%e, %e, %e)' )';
      end
    end

    data(N_data).variables = ln;
    num_col = 14;
    if strcmp(data(N_data).variables,'t Tp1 Tp2 Tp3 Tn V Vp1 Vp2 Vp3 Dp Dp1 Dp2 Dp3 Dn') ~= 1
      error('Format changed')
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

  pdata.t = data(1).t(1:n_steps);


  pdata.Dp  = 0;
  pdata.Dp1 = 0;
  pdata.Dp2 = 0;
  pdata.Dp3 = 0;
  pdata.Dn  = 0;

  pdata.V   = 0;
  pdata.Vp1 = 0;
  pdata.Vp2 = 0;
  pdata.Vp3 = 0;

  pdata.Tn  = 0;
  pdata.Tp1 = 0;
  pdata.Tp2 = 0;
  pdata.Tp3 = 0;

  pdata.n   = sign(n_ref(1)*data(1).n(1) + n_ref(2)*data(1).n(2) + n_ref(3)*data(1).n(3))*data(1).n;
  head_saved = false;
  for n = 1:N_data
    n_sign    = sign(n_ref(1)*data(n).n(1) + n_ref(2)*data(n).n(2) + n_ref(3)*data(n).n(3));
    if(~head_saved && n_sign == 1)
      head_saved = true;
      pdata.header = data(n).header;
    end
    if(pdata.n ~= n_sign*data(n).n)
      warning(['Normals maybe discontinuous (rethink averaging)!'...
               'May want to work with tractions directly????'])
      disp(pdata.n)
      disp(n_sign*data(n).n)
    end

    pdata.Dp  = pdata.Dp  +        data(n).Dp /N_data;
    disp(size(pdata.Dp1))
    disp(size(data(n).Dp1))
    pdata.Dp1 = pdata.Dp1 + n_sign*data(n).Dp1/N_data;
    pdata.Dp2 = pdata.Dp2 + n_sign*data(n).Dp2/N_data;
    pdata.Dp3 = pdata.Dp3 + n_sign*data(n).Dp3/N_data;

    pdata.V   = pdata.V   +        data(n).V  /N_data;
    pdata.Vp1 = pdata.Vp1 + n_sign*data(n).Vp1/N_data;
    pdata.Vp2 = pdata.Vp2 + n_sign*data(n).Vp2/N_data;
    pdata.Vp3 = pdata.Vp3 + n_sign*data(n).Vp3/N_data;

    pdata.Tn  = pdata.Tn  +        data(n).Tn /N_data;
    pdata.Tp1 = pdata.Tp1 + n_sign*data(n).Tp1/N_data;
    pdata.Tp2 = pdata.Tp2 + n_sign*data(n).Tp2/N_data;
    pdata.Tp3 = pdata.Tp3 + n_sign*data(n).Tp3/N_data;

  end
  pdata.Tn (1) = pdata.Tn (2) ;
  pdata.Tp1(1) = pdata.Tp1(2);
  pdata.Tp2(1) = pdata.Tp2(2);
  pdata.Tp3(1) = pdata.Tp3(2);
else
  pdata = [];
end
