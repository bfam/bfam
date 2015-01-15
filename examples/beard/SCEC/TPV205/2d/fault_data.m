clear

fault_stations = { 'm4.5_0', 'm7.5_0', 'm12_0' , '0_0'   , 'p4.5_0', 'p7.5_0', 'p12_0'};

for k = 1:length(fault_stations)
  disp(fault_stations{k})
  base = ['TPV205_rot_pi_4_',fault_stations{k}];
  data = process_fault_station_2d(base, ['scec/',base,'.scec'], 'data');
end
