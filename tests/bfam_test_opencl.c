#include <bfam.h>
#include <bfam_test_opencl_kernels.h>

#define SIZE 34

int
main(int argc, char *argv[])
{
  int             failures = 0;
#ifdef BFAM_HAVE_OPENCL

  cl_context      context;
  cl_int          status;
  cl_command_queue queue;

  bfam_cl_print_platforms_devices();
  bfam_cl_create_context_on(NULL, NULL, 0, &context, &queue, 0);
  cl_kernel kernel = bfam_cl_kernel_from_string(context,
                                                bfam_test_opencl_kernels,
                                                "square", NULL);

  bfam_cl_print_device_info_from_queue(queue);

  cl_int * h_data = bfam_calloc(SIZE, sizeof(cl_int));
  for (size_t i = 0; i < SIZE; ++i) {
    h_data[i] = i;
  }

  /*
   * Allocate memory
   */
  cl_mem c_data = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR,
                                 sizeof(cl_int) * SIZE, h_data, &status);
  BFAM_CL_CHECK(status, "clCreateBuffer");

  /*
   * Set the argument to the kernel
   */
  BFAM_CL_SET_1_KERNEL_ARG(kernel, c_data);

  size_t ldim[] = {SIZE};
  size_t gdim[] = {SIZE, 0, 0};

  /*
   * submit kernel to the queue
   */
  BFAM_CL_SAFE_CALL(
      clEnqueueNDRangeKernel(queue, kernel, 1, NULL, gdim,
                             ldim, 0, NULL, NULL)
      );

  /*
   * transfer the result to the host
   */
  BFAM_CL_SAFE_CALL(
      clEnqueueReadBuffer(queue, c_data, CL_TRUE, 0,
                          sizeof(cl_int) * SIZE,
                          h_data, 0, NULL, NULL)
      );

  /*
   * verify the result
   */
  for (size_t i = 0; i < SIZE; ++i) {
    BFAM_INFO("%d: %d", (int) i, (int) h_data[i]);
    if (h_data[i] != (cl_int) (i * i))
      ++failures;
  }

  /*
   * clean up
   */
  BFAM_CL_SAFE_CALL(clReleaseMemObject(c_data));
  BFAM_CL_SAFE_CALL(clReleaseKernel(kernel));
  BFAM_CL_SAFE_CALL(clReleaseCommandQueue(queue));
  BFAM_CL_SAFE_CALL(clReleaseContext(context));
  bfam_free(h_data);

#else
  BFAM_WARNING("BFAM not compiled with OpenCL support");
#endif

  return failures;
}
