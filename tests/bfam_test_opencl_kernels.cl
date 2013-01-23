void kernel square(global int *data) {
  int id = get_global_id(0);
  data[id] *= data[id];
}
