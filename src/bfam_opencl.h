/*
 * The functions in this file are modified versions of the ones found in
 * Andreas Kloeckner's cl-helper routines from:
 *
 *   https://github.com/hpc12/lec1-demo.git
 *
 * See the copyright statement below for license information.
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

#ifndef BFAM_OPENCL_H
#define BFAM_OPENCL_H

#include <bfam_base.h>
#include <bfam_log.h>

#ifdef BFAM_HAVE_OPENCL

#define BFAM_CL_CHECK(err, call) do                                            \
  {                                                                            \
    if ((err) != CL_SUCCESS )                                                  \
    {                                                                          \
      BFAM_LERROR("OpenCL Error: %s [%s]", bfam_cl_err_to_str((err)), (call)); \
      BFAM_ABORT("OpenCL Error");                                              \
    }                                                                          \
  } while (0)

#define BFAM_CL_SAFE_CALL(call) do                  \
  {                                                 \
    cl_int bfam_cl_err_xxxuniquexxx = call;         \
    BFAM_CL_CHECK(bfam_cl_err_xxxuniquexxx, #call); \
  } while (0)

/** Get error string associated with an OpenCL error code.
 *
 * \param[in] err OpenCL error code
 * \return A pointer to a string literal containing the error string.
 */
const char * bfam_cl_err_to_str(cl_int e);

/** Print OpenCL Platforms and Devices.
 */
void bfam_cl_print_platforms_devices();

/** Create an OpenCL context and a matching command queue.
 *
 * Create an OpenCL context and a matching command queue on a platform from a
 * vendor whose name contains \a plat_name on a device whose name contains
 * \a dev_name. Both \a plat_name and \a dev_name may be \c NULL, indicating no
 * preference in the matter.
 *
 * If multiple devices match both \a plat_name and \a dev_name, then \a idx
 * prescribes the number of the device that should be chosen.
 *
 * This function always succeeds. (If an error occurs, the program
 * is aborted.)
 *
 * \param[in]  plat_name        platform name (\c NULL indicates no preference)
 * \param[in]  dev_name         device name (\c NULL indicates no preference)
 * \param[in]  idx              the number of the device that should be chosen
 * \param[out] ctx              returned OpenCL context
 * \param[out] queue            returned OpenCL queue
 * \param[in]  enable_profiling boolean indicating to turn on profiling
 */
void bfam_cl_create_context_on(const char *plat_name, const char *dev_name,
                               cl_uint idx, cl_context *ctx,
                               cl_command_queue *queue, int enable_profiling);


/** Create an OpenCL kernel.
 *
 * Create a new OpenCL kernel from the code in the string \a knl.
 * \a knl_name is the name of the kernel function, and \a options,
 * if not \c NULL, is a string containing compiler flags.
 *
 * You must release the resulting kernel when you're done
 * with it.
 *
 * This function always succeeds. (If an error occurs, the program
 * is aborted.)
 *
 * \param[in] ctx
 * \param[in] knl
 * \param[in] knl_name
 * \param[in] options
 * \returns Open
 */
cl_kernel bfam_cl_kernel_from_string(cl_context ctx, char const *knl,
                                     char const *knl_name, char const *options);

/** Print information about a device.
 *
 * \param[in] device OpenCL device to print information about
 */
void bfam_cl_print_device_info(cl_device_id device);

/** Print information about a device from a queue.
 *
 * \param[in] queue OpenCL queue to print information about the associated
 *                  device
 */
void bfam_cl_print_device_info_from_queue(cl_command_queue queue);

#define BFAM_CL_SET_1_KERNEL_ARG(knl, arg0) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0));

#define BFAM_CL_SET_2_KERNEL_ARGS(knl, arg0, arg1) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1));

#define BFAM_CL_SET_3_KERNEL_ARGS(knl, arg0, arg1, arg2) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2));

#define BFAM_CL_SET_4_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3));

#define BFAM_CL_SET_5_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4));

#define BFAM_CL_SET_6_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5));

#define BFAM_CL_SET_7_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5, arg6) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 6, sizeof(arg6), &arg6));

#define BFAM_CL_SET_8_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 6, sizeof(arg6), &arg6)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 7, sizeof(arg7), &arg7));

#define BFAM_CL_SET_9_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 6, sizeof(arg6), &arg6)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 7, sizeof(arg7), &arg7)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 8, sizeof(arg8), &arg8));

#define BFAM_CL_SET_10_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 6, sizeof(arg6), &arg6)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 7, sizeof(arg7), &arg7)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 8, sizeof(arg8), &arg8)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 9, sizeof(arg9), &arg9));

#define BFAM_CL_SET_11_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 6, sizeof(arg6), &arg6)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 7, sizeof(arg7), &arg7)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 8, sizeof(arg8), &arg8)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 9, sizeof(arg9), &arg9)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 10, sizeof(arg10), &arg10));

#define BFAM_CL_SET_12_KERNEL_ARGS(knl, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11) \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 0, sizeof(arg0), &arg0)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 1, sizeof(arg1), &arg1)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 2, sizeof(arg2), &arg2)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 3, sizeof(arg3), &arg3)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 4, sizeof(arg4), &arg4)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 5, sizeof(arg5), &arg5)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 6, sizeof(arg6), &arg6)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 7, sizeof(arg7), &arg7)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 8, sizeof(arg8), &arg8)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 9, sizeof(arg9), &arg9)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 10, sizeof(arg10), &arg10)); \
  BFAM_CL_SAFE_CALL(clSetKernelArg(knl, 11, sizeof(arg11), &arg11));



#endif

#endif
