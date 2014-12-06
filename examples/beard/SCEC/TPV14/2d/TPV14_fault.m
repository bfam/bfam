clear

sim_name = 'TPV14_base_';
names = {[sim_name,'faultst-020dp075'],...
         [sim_name,'faultst020dp075'] ,...
         [sim_name,'faultst050dp075'] ,...
         [sim_name,'faultst090dp075'] ,...
         [sim_name,'branchst020dp075'],...
         [sim_name,'branchst050dp075'],...
         [sim_name,'branchst090dp075']};

for k = 1:length(names)
  pd = process_fault_station([0,1],names{k},'data');

  if(isempty(pd))
    disp(['NOT FOUND :: ',names{k}])
  else
    disp(['    FOUND :: ',names{k}])
    m = [pd.n(2);-pd.n(1)];
    t       = pd.t;
    hslip   = pd.Dp1*m(1)+pd.Dp2*m(2);
    hrate   = pd.Vp1*m(1)+pd.Vp2*m(2);
    hshear  = pd.Tp1*m(1)+pd.Tp2*m(2);
    vslip   = zeros(size(hslip));
    vrate   = zeros(size(hrate));
    vshear  = zeros(size(hshear));
    nstress = pd.Tn;
    write_scec_data(['scec/',names{k},'.scec'],pd.header,...
       't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress',...
       [t, hslip, hrate, hshear, vslip, vrate, vshear, nstress]);
  end
end
