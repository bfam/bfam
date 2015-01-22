#include <bfam_base.h>
#include <bfam_opencl.h>

/*
 * The functions in this file are modified versions of the ones found in
 * Andreas Kloeckner's cl-helper routines.  See the copyright statement below
 * for license information.
 *
 * Copyright (c) 2010, 2012 Andreas Kloeckner
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifdef BFAM_HAVE_OPENCL

#define MAX_NAME_LEN BUFSIZ

const char *bfam_cl_err_to_str(cl_int e)
{
  switch (e)
  {
  case CL_SUCCESS:
    return "success";
  case CL_DEVICE_NOT_FOUND:
    return "device not found";
  case CL_DEVICE_NOT_AVAILABLE:
    return "device not available";
#if !(defined(CL_PLATFORM_NVIDIA) && CL_PLATFORM_NVIDIA == 0x3001)
  case CL_COMPILER_NOT_AVAILABLE:
    return "device compiler not available";
#endif
  case CL_MEM_OBJECT_ALLOCATION_FAILURE:
    return "mem object allocation failure";
  case CL_OUT_OF_RESOURCES:
    return "out of resources";
  case CL_OUT_OF_HOST_MEMORY:
    return "out of host memory";
  case CL_PROFILING_INFO_NOT_AVAILABLE:
    return "profiling info not available";
  case CL_MEM_COPY_OVERLAP:
    return "mem copy overlap";
  case CL_IMAGE_FORMAT_MISMATCH:
    return "image format mismatch";
  case CL_IMAGE_FORMAT_NOT_SUPPORTED:
    return "image format not supported";
  case CL_BUILD_PROGRAM_FAILURE:
    return "build program failure";
  case CL_MAP_FAILURE:
    return "map failure";

  case CL_INVALID_VALUE:
    return "invalid value";
  case CL_INVALID_DEVICE_TYPE:
    return "invalid device type";
  case CL_INVALID_PLATFORM:
    return "invalid platform";
  case CL_INVALID_DEVICE:
    return "invalid device";
  case CL_INVALID_CONTEXT:
    return "invalid context";
  case CL_INVALID_QUEUE_PROPERTIES:
    return "invalid queue properties";
  case CL_INVALID_COMMAND_QUEUE:
    return "invalid command queue";
  case CL_INVALID_HOST_PTR:
    return "invalid host ptr";
  case CL_INVALID_MEM_OBJECT:
    return "invalid mem object";
  case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
    return "invalid image format descriptor";
  case CL_INVALID_IMAGE_SIZE:
    return "invalid image size";
  case CL_INVALID_SAMPLER:
    return "invalid sampler";
  case CL_INVALID_BINARY:
    return "invalid binary";
  case CL_INVALID_BUILD_OPTIONS:
    return "invalid build options";
  case CL_INVALID_PROGRAM:
    return "invalid program";
  case CL_INVALID_PROGRAM_EXECUTABLE:
    return "invalid program executable";
  case CL_INVALID_KERNEL_NAME:
    return "invalid kernel name";
  case CL_INVALID_KERNEL_DEFINITION:
    return "invalid kernel definition";
  case CL_INVALID_KERNEL:
    return "invalid kernel";
  case CL_INVALID_ARG_INDEX:
    return "invalid arg index";
  case CL_INVALID_ARG_VALUE:
    return "invalid arg value";
  case CL_INVALID_ARG_SIZE:
    return "invalid arg size";
  case CL_INVALID_KERNEL_ARGS:
    return "invalid kernel args";
  case CL_INVALID_WORK_DIMENSION:
    return "invalid work dimension";
  case CL_INVALID_WORK_GROUP_SIZE:
    return "invalid work group size";
  case CL_INVALID_WORK_ITEM_SIZE:
    return "invalid work item size";
  case CL_INVALID_GLOBAL_OFFSET:
    return "invalid global offset";
  case CL_INVALID_EVENT_WAIT_LIST:
    return "invalid event wait list";
  case CL_INVALID_EVENT:
    return "invalid event";
  case CL_INVALID_OPERATION:
    return "invalid operation";
  case CL_INVALID_GL_OBJECT:
    return "invalid gl object";
  case CL_INVALID_BUFFER_SIZE:
    return "invalid buffer size";
  case CL_INVALID_MIP_LEVEL:
    return "invalid mip level";

#if defined(cl_khr_gl_sharing) && (cl_khr_gl_sharing >= 1)
  case CL_INVALID_GL_SHAREGROUP_REFERENCE_KHR:
    return "invalid gl sharegroup reference number";
#endif

#ifdef CL_VERSION_1_1
  case CL_MISALIGNED_SUB_BUFFER_OFFSET:
    return "misaligned sub-buffer offset";
  case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
    return "exec status error for events in wait list";
  case CL_INVALID_GLOBAL_WORK_SIZE:
    return "invalid global work size";
#endif

  default:
    return "invalid/unknown error code";
  }
}

void bfam_cl_print_platforms_devices()
{
  // get number of platforms
  cl_uint plat_count;
  BFAM_CL_SAFE_CALL(clGetPlatformIDs(0, NULL, &plat_count));

  // allocate memory, get list of platforms
  cl_platform_id *platforms =
      (cl_platform_id *)bfam_malloc(plat_count * sizeof(cl_platform_id));

  BFAM_CL_SAFE_CALL(clGetPlatformIDs(plat_count, platforms, NULL));

  for (cl_uint i = 0; i < plat_count; ++i)
  {
    // get platform vendor name
    char buf[MAX_NAME_LEN];
    BFAM_CL_SAFE_CALL(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
                                        sizeof(buf), buf, NULL));
    BFAM_INFO("platform %d: vendor '%s'", i, buf);

    // get number of devices in platform
    cl_uint dev_count;
    BFAM_CL_SAFE_CALL(
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &dev_count));

    cl_device_id *devices =
        (cl_device_id *)bfam_malloc(dev_count * sizeof(cl_device_id));

    // get list of devices in platform
    BFAM_CL_SAFE_CALL(clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL,
                                     dev_count, devices, NULL));

    // iterate over devices
    for (cl_uint j = 0; j < dev_count; ++j)
    {
      char buf[MAX_NAME_LEN];
      BFAM_CL_SAFE_CALL(
          clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(buf), buf, NULL));
      BFAM_INFO("  device %d: '%s'", j, buf);
    }

    bfam_free(devices);
  }

  bfam_free(platforms);
}

/* Read a line from stdin. C makes things simple. :)
 * From http://stackoverflow.com/a/314422/1148634
 */
static char *read_a_line(void)
{
  char *line = (char *)bfam_malloc(MAX_NAME_LEN), *linep = line;
  size_t lenmax = MAX_NAME_LEN, len = lenmax;
  int c;

  if (line == NULL)
    return NULL;

  for (;;)
  {
    c = fgetc(stdin);
    if (c == EOF)
      break;

    if (--len == 0)
    {
      char *linen = (char *)bfam_realloc(linep, lenmax *= 2);
      len = lenmax;

      if (linen == NULL)
      {
        bfam_free(linep);
        return NULL;
      }
      line = linen + (line - linep);
      linep = linen;
    }

    if ((*line++ = c) == '\n')
      break;
  }
  *line = '\0';
  return linep;
}

const char *CHOOSE_INTERACTIVELY = "INTERACTIVE";

void bfam_cl_create_context_on(const char *plat_name, const char *dev_name,
                               cl_uint idx, cl_context *ctx,
                               cl_command_queue *queue, int enable_profiling)
{
  char dev_sel_buf[MAX_NAME_LEN];
  char platform_sel_buf[MAX_NAME_LEN];

  // get number of platforms
  cl_uint plat_count;
  BFAM_CL_SAFE_CALL(clGetPlatformIDs(0, NULL, &plat_count));

  // allocate memory, get list of platform handles
  cl_platform_id *platforms =
      (cl_platform_id *)bfam_malloc(plat_count * sizeof(cl_platform_id));
  BFAM_CL_SAFE_CALL(clGetPlatformIDs(plat_count, platforms, NULL));

// print menu, if requested
#ifndef BFAM_CL_HELPER_FORCE_INTERACTIVE
  if (plat_name == CHOOSE_INTERACTIVELY) // yes, we want exactly that pointer
#endif
  {
    puts("Choose platform:");
    for (cl_uint i = 0; i < plat_count; ++i)
    {
      char buf[MAX_NAME_LEN];
      BFAM_CL_SAFE_CALL(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
                                          sizeof(buf), buf, NULL));
      printf("[%d] %s\n", i, buf);
    }

    printf("Enter choice: ");
    fflush(stdout);

    char *sel = read_a_line();
    if (!sel)
    {
      BFAM_ABORT("error reading line from stdin");
    }

    int sel_int = BFAM_MIN(BFAM_MAX(0, atoi(sel)), (int)plat_count - 1);
    bfam_free(sel);

    BFAM_CL_SAFE_CALL(clGetPlatformInfo(platforms[sel_int], CL_PLATFORM_VENDOR,
                                        sizeof(platform_sel_buf),
                                        platform_sel_buf, NULL));
    plat_name = platform_sel_buf;
  }

  // iterate over platforms
  for (cl_uint i = 0; i < plat_count; ++i)
  {
    // get platform name
    char buf[MAX_NAME_LEN];
    BFAM_CL_SAFE_CALL(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
                                        sizeof(buf), buf, NULL));

    // does it match?
    if (!plat_name || strstr(buf, plat_name))
    {
      // get number of devices in platform
      cl_uint dev_count;
      BFAM_CL_SAFE_CALL(clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0,
                                       NULL, &dev_count));

      // allocate memory, get list of device handles in platform
      cl_device_id *devices =
          (cl_device_id *)bfam_malloc(dev_count * sizeof(cl_device_id));

      BFAM_CL_SAFE_CALL(clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL,
                                       dev_count, devices, NULL));

// {{{ print device menu, if requested
#ifndef BFAM_CL_HELPER_FORCE_INTERACTIVE
      if (dev_name == CHOOSE_INTERACTIVELY) // yes, we want exactly that pointer
#endif
      {
        puts("Choose device:");
        for (cl_uint j = 0; j < dev_count; ++j)
        {
          char buf[MAX_NAME_LEN];
          BFAM_CL_SAFE_CALL(clGetDeviceInfo(devices[j], CL_DEVICE_NAME,
                                            sizeof(buf), buf, NULL));
          printf("[%d] %s\n", j, buf);
        }

        printf("Enter choice: ");
        fflush(stdout);

        char *sel = read_a_line();
        if (!sel)
        {
          BFAM_ABORT("error reading line from stdin");
          abort(); /* Add this for the clang static analysis */
        }

        int int_sel = BFAM_MIN(BFAM_MAX(0, atoi(sel)), (int)dev_count - 1);
        bfam_free(sel);

        BFAM_CL_SAFE_CALL(clGetDeviceInfo(devices[int_sel], CL_DEVICE_NAME,
                                          sizeof(dev_sel_buf), dev_sel_buf,
                                          NULL));
        dev_name = dev_sel_buf;
      }

      // }}}

      // iterate over devices
      for (cl_uint j = 0; j < dev_count; ++j)
      {
        // get device name
        char buf[MAX_NAME_LEN];
        BFAM_CL_SAFE_CALL(clGetDeviceInfo(devices[j], CL_DEVICE_NAME,
                                          sizeof(buf), buf, NULL));

        // does it match?
        if (!dev_name || strstr(buf, dev_name))
        {
          if (idx == 0)
          {
            cl_platform_id plat = platforms[i];
            cl_device_id dev = devices[j];

            bfam_free(devices);
            bfam_free(platforms);

            // create a context
            cl_context_properties cps[3] = {CL_CONTEXT_PLATFORM,
                                            (cl_context_properties)plat, 0};

            cl_int status;
            *ctx = clCreateContext(cps, 1, &dev, NULL, NULL, &status);
            BFAM_CL_CHECK(status, "clCreateContext");

            // create a command queue
            cl_command_queue_properties qprops = 0;
            if (enable_profiling)
              qprops |= CL_QUEUE_PROFILING_ENABLE;

            *queue = clCreateCommandQueue(*ctx, dev, qprops, &status);
            BFAM_CL_CHECK(status, "clCreateCommandQueue");

            return;
          }
          else
            --idx;
        }
      }

      bfam_free(devices);
    }
  }

  bfam_free(platforms);

  BFAM_ABORT("create_context_on: specified device not found.");
}

cl_kernel bfam_cl_kernel_from_string(cl_context ctx, char const *knl,
                                     char const *knl_name, char const *options)
{
  // create an OpenCL program (may have multiple kernels)
  size_t sizes[] = {strlen(knl)};

  cl_int status;
  cl_program program = clCreateProgramWithSource(ctx, 1, &knl, sizes, &status);
  BFAM_CL_CHECK(status, "clCreateProgramWithSource");

  // build it
  status = clBuildProgram(program, 0, NULL, options, NULL, NULL);

  if (status != CL_SUCCESS)
  {
    // build failed, get build log and print it

    cl_device_id dev;
    BFAM_CL_SAFE_CALL(
        clGetProgramInfo(program, CL_PROGRAM_DEVICES, sizeof(dev), &dev, NULL));

    size_t log_size;
    BFAM_CL_SAFE_CALL(clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
                                            0, NULL, &log_size));

    char *log = (char *)bfam_malloc(log_size);

    char devname[MAX_NAME_LEN];
    BFAM_CL_SAFE_CALL(
        clGetDeviceInfo(dev, CL_DEVICE_NAME, sizeof(devname), devname, NULL));

    BFAM_CL_SAFE_CALL(clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
                                            log_size, log, NULL));
    BFAM_LERROR("*** build of '%s' on '%s' failed:\n%s\n*** (end of error)\n",
                knl_name, devname, log);
    BFAM_ABORT("Building kernel from a string");
  }
  else
    BFAM_CL_CHECK(status, "clBuildProgram");

  // fish the kernel out of the program
  cl_kernel kernel = clCreateKernel(program, knl_name, &status);
  BFAM_CL_CHECK(status, "clCreateKernel");

  BFAM_CL_SAFE_CALL(clReleaseProgram(program));

  return kernel;
}

void bfam_cl_print_device_info(cl_device_id device)
{
// adapted from http://graphics.stanford.edu/~yoel/notes/clInfo.c

#define LONG_PROPS                                                             \
  defn(VENDOR_ID), defn(MAX_COMPUTE_UNITS), defn(MAX_WORK_ITEM_DIMENSIONS),    \
      defn(MAX_WORK_GROUP_SIZE), defn(PREFERRED_VECTOR_WIDTH_CHAR),            \
      defn(PREFERRED_VECTOR_WIDTH_SHORT), defn(PREFERRED_VECTOR_WIDTH_INT),    \
      defn(PREFERRED_VECTOR_WIDTH_LONG), defn(PREFERRED_VECTOR_WIDTH_FLOAT),   \
      defn(PREFERRED_VECTOR_WIDTH_DOUBLE), defn(MAX_CLOCK_FREQUENCY),          \
      defn(ADDRESS_BITS), defn(MAX_MEM_ALLOC_SIZE), defn(IMAGE_SUPPORT),       \
      defn(MAX_READ_IMAGE_ARGS), defn(MAX_WRITE_IMAGE_ARGS),                   \
      defn(IMAGE2D_MAX_WIDTH), defn(IMAGE2D_MAX_HEIGHT),                       \
      defn(IMAGE3D_MAX_WIDTH), defn(IMAGE3D_MAX_HEIGHT),                       \
      defn(IMAGE3D_MAX_DEPTH), defn(MAX_SAMPLERS), defn(MAX_PARAMETER_SIZE),   \
      defn(MEM_BASE_ADDR_ALIGN), defn(MIN_DATA_TYPE_ALIGN_SIZE),               \
      defn(GLOBAL_MEM_CACHELINE_SIZE), defn(GLOBAL_MEM_CACHE_SIZE),            \
      defn(GLOBAL_MEM_SIZE), defn(MAX_CONSTANT_BUFFER_SIZE),                   \
      defn(MAX_CONSTANT_ARGS), defn(LOCAL_MEM_SIZE),                           \
      defn(ERROR_CORRECTION_SUPPORT), defn(PROFILING_TIMER_RESOLUTION),        \
      defn(ENDIAN_LITTLE), defn(AVAILABLE), defn(COMPILER_AVAILABLE),

#define STR_PROPS                                                              \
  defn(NAME), defn(VENDOR), defn(PROFILE), defn(VERSION), defn(EXTENSIONS),

#define HEX_PROPS defn(SINGLE_FP_CONFIG), defn(QUEUE_PROPERTIES),

  printf("---------------------------------------------------------------------"
         "\n");

  static struct
  {
    cl_device_info param;
    const char *name;
  } longProps[] = {
#define defn(X)                                                                \
  {                                                                            \
    CL_DEVICE_##X, #X                                                          \
  }
        LONG_PROPS
#undef defn
        {0, NULL},
    };
  static struct
  {
    cl_device_info param;
    const char *name;
  } hexProps[] = {
#define defn(X)                                                                \
  {                                                                            \
    CL_DEVICE_##X, #X                                                          \
  }
        HEX_PROPS
#undef defn
        {0, NULL},
    };
  static struct
  {
    cl_device_info param;
    const char *name;
  } strProps[] = {
#define defn(X)                                                                \
  {                                                                            \
    CL_DEVICE_##X, #X                                                          \
  }
        STR_PROPS
#undef defn
        {CL_DRIVER_VERSION, "DRIVER_VERSION"},
        {0, NULL},
    };
  cl_int status;
  size_t size;
  char buf[65536];
  long long val; /* Avoids unpleasant surprises for some params */
  int ii;

  for (ii = 0; strProps[ii].name != NULL; ii++)
  {
    status =
        clGetDeviceInfo(device, strProps[ii].param, sizeof buf, buf, &size);
    if (status != CL_SUCCESS)
    {
      printf("Unable to get %s: %s!\n", strProps[ii].name,
             bfam_cl_err_to_str(status));
      continue;
    }
    if (size > sizeof buf)
    {
      printf("Large %s (%zd bytes)!  Truncating to %ld!\n", strProps[ii].name,
             size, sizeof buf);
    }
    printf("%s: %s\n", strProps[ii].name, buf);
  }
  printf("\n");

  status = clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof val, &val, NULL);
  if (status == CL_SUCCESS)
  {
    printf("Type: ");
    if (val & CL_DEVICE_TYPE_DEFAULT)
    {
      val &= ~CL_DEVICE_TYPE_DEFAULT;
      printf("Default ");
    }
    if (val & CL_DEVICE_TYPE_CPU)
    {
      val &= ~CL_DEVICE_TYPE_CPU;
      printf("CPU ");
    }
    if (val & CL_DEVICE_TYPE_GPU)
    {
      val &= ~CL_DEVICE_TYPE_GPU;
      printf("GPU ");
    }
    if (val & CL_DEVICE_TYPE_ACCELERATOR)
    {
      val &= ~CL_DEVICE_TYPE_ACCELERATOR;
      printf("Accelerator ");
    }
    if (val != 0)
    {
      printf("Unknown (0x%llx) ", val);
    }
    printf("\n");
  }
  else
  {
    printf("Unable to get TYPE: %s!\n", bfam_cl_err_to_str(status));
  }

  status = clGetDeviceInfo(device, CL_DEVICE_EXECUTION_CAPABILITIES, sizeof val,
                           &val, NULL);
  if (status == CL_SUCCESS)
  {
    printf("EXECUTION_CAPABILITIES: ");
    if (val & CL_EXEC_KERNEL)
    {
      val &= ~CL_EXEC_KERNEL;
      printf("Kernel ");
    }
    if (val & CL_EXEC_NATIVE_KERNEL)
    {
      val &= ~CL_EXEC_NATIVE_KERNEL;
      printf("Native ");
    }
    if (val)
      printf("Unknown (0x%llx) ", val);

    printf("\n");
  }
  else
  {
    printf("Unable to get EXECUTION_CAPABILITIES: %s!\n",
           bfam_cl_err_to_str(status));
  }

  status = clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHE_TYPE, sizeof val,
                           &val, NULL);
  if (status == CL_SUCCESS)
  {
    static const char *cacheTypes[] = {"None", "Read-Only", "Read-Write"};
    static int numTypes = sizeof cacheTypes / sizeof cacheTypes[0];

    printf("GLOBAL_MEM_CACHE_TYPE: %s (%lld)\n",
           val < numTypes ? cacheTypes[val] : "???", val);
  }
  else
  {
    printf("Unable to get GLOBAL_MEM_CACHE_TYPE: %s!\n",
           bfam_cl_err_to_str(status));
  }

  status =
      clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof val, &val, NULL);

  if (status == CL_SUCCESS)
  {
    static const char *lmemTypes[] = {"???", "Local", "Global"};
    static int numTypes = sizeof lmemTypes / sizeof lmemTypes[0];

    printf("CL_DEVICE_LOCAL_MEM_TYPE: %s (%lld)\n",
           val < numTypes ? lmemTypes[val] : "???", val);
  }
  else
  {
    printf("Unable to get CL_DEVICE_LOCAL_MEM_TYPE: %s!\n",
           bfam_cl_err_to_str(status));
  }

  for (ii = 0; hexProps[ii].name != NULL; ii++)
  {
    status =
        clGetDeviceInfo(device, hexProps[ii].param, sizeof val, &val, &size);
    if (status != CL_SUCCESS)
    {
      printf("Unable to get %s: %s!\n", hexProps[ii].name,
             bfam_cl_err_to_str(status));
      continue;
    }
    if (size > sizeof val)
    {
      printf("Large %s (%zd bytes)!  Truncating to %ld!\n", hexProps[ii].name,
             size, sizeof val);
    }
    printf("%s: 0x%llx\n", hexProps[ii].name, val);
  }
  printf("\n");

  for (ii = 0; longProps[ii].name != NULL; ii++)
  {
    status =
        clGetDeviceInfo(device, longProps[ii].param, sizeof val, &val, &size);
    if (status != CL_SUCCESS)
    {
      printf("Unable to get %s: %s!\n", longProps[ii].name,
             bfam_cl_err_to_str(status));
      continue;
    }
    if (size > sizeof val)
    {
      printf("Large %s (%zd bytes)!  Truncating to %ld!\n", longProps[ii].name,
             size, sizeof val);
    }
    printf("%s: %lld\n", longProps[ii].name, val);
  }

  {
    size_t size;
    BFAM_CL_SAFE_CALL(
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, 0, 0, &size));

    size_t res_vec[size / sizeof(size_t)]; // C99 VLA yay!

    BFAM_CL_SAFE_CALL(clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES,
                                      size, res_vec, &size));

    printf("MAX_WORK_GROUP_SIZES: "); // a tiny lie
    for (size_t i = 0; i < size / sizeof(size_t); ++i)
      printf("%zd ", res_vec[i]);
    printf("\n");
  }
  printf("---------------------------------------------------------------------"
         "\n");
}

void bfam_cl_print_device_info_from_queue(cl_command_queue queue)
{
  cl_device_id dev;
  BFAM_CL_SAFE_CALL(
      clGetCommandQueueInfo(queue, CL_QUEUE_DEVICE, sizeof dev, &dev, NULL));

  bfam_cl_print_device_info(dev);
}

#endif
