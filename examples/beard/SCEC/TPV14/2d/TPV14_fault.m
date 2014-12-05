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
    plot(pd.t,pd.V)
    pause
  end
  write_scec_data('tmp',pd.header)
end
