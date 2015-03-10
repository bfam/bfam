# Patch for OpenCL on Yosemite

To get OpenCL working with gcc you need to apply a patch to the system headers.
This can be done with the command:

```
sudo patch -b -N -p0 <dispatch.patch
```

This patch came from this thread:

  <https://github.com/magnumripper/JohnTheRipper/issues/1025>
