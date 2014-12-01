function data = process_fault_station(base_name, outputname, directory)

if nargin < 3
  directory = '.';
end
if nargin < 2
  outputname = [base_name,'.scec'];
end

A = dir(directory);

data = [];

cnt = 0;
for k = 1:length(A)
  if strncmpi(base_name,A(k).name,length(base_name)) && 1 ~=  strcmp(outputname,A(k).name)
    cnt = cnt + 1;
    data(cnt).name = A(k).name;

    full_name = [directory,'/',A(k).name];

    % Find header lines
    iid = fopen(full_name, 'r');
    header = 0;
    ln = fgetl(iid);
    while ln(1) == '#'
      header = header + 1;
      ln = fgetl(iid);
      if strcmp(ln,'# interpolated normal:') == 1
        header = header + 1;
        ln = fgetl(iid);
        data(cnt).n = sscanf(ln,'#    (%e, %e, %e)' )';
      end
    end

    data(cnt).variables = ln;
    num_col = 14;
    if strcmp(data(cnt).variables,'t Tp1 Tp2 Tp3 Tn V Vp1 Vp2 Vp3 Dp Dp1 Dp2 Dp3 Dn') ~= 1
      error('Format changed')
    end
    dat = fscanf(iid,'%e');
    data(cnt).dat = reshape(dat,num_col,length(dat)/num_col)';

    fclose(iid);
  end
end

[t,Tp1,Tp2,Tp3,Tn,V,Vp1,Vp2,Vp3,Dp,Dp1,Dp2,Dp3,Dn] = deal(1,2,3,4,5,6,7,8,9,10,11,12,13,14);
if cnt > 0
  M = length(data(1).dat(:,t));
  for k = 2:cnt
    M = min(M,length(data(k).dat(:,t)));
  end
  t = data(1).dat(1:M,t);
  hslip        = zeros(size(t));
  hsliprate    = zeros(size(t));
  hshearstress = zeros(size(t));
  vslip        = zeros(size(t));
  vsliprate    = zeros(size(t));
  vshearstress = zeros(size(t));
  nstress      = zeros(size(t));
  for k = 1:cnt
    vshearstress = vshearstress - data(cnt).dat(1:M,Tp3)/cnt/data(cnt).n(2);
    vsliprate    = vsliprate    + data(cnt).dat(1:M,Vp3)/cnt;
    vslip        = vslip        - data(cnt).dat(1:M,Dp3)/cnt;

    hshearstress = hshearstress + data(cnt).dat(1:M,Tp1)/cnt/data(cnt).n(2);;
    hsliprate    = hsliprate    + data(cnt).dat(1:M,Vp1)/cnt;
    hslip        = hslip        + data(cnt).dat(1:M,Dp1)/cnt;

    nstress      = nstress      + data(cnt).dat(1:M,Tn )/cnt;
  end

  oid = fopen(outputname,'w');
  iid = fopen(full_name, 'r');
  ln = fgetl(iid);
  while ln(1) == '#'
    fprintf(oid,'%s\n',ln);
    ln = fgetl(iid);
  end
  fclose(iid);

  fprintf(oid,'t  h-slip  h-slip-rate  h-shear-stress  v-slip  v-slip-rate v-shear-stress\n');
  for k = 2:length(t)
    fprintf(oid,'%+30.17E %+30.17E %+30.17E %+30.17E %+30.17E %+30.17E %+30.17E\n',...
             t(k),hslip(k),hsliprate(k),hshearstress(k),vslip(k),vsliprate(k),vshearstress(k));
  end
  fclose(oid);
end


