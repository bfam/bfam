/*
 * This code is adapted from:
 *
 *   https://github.com/ThomasHabets/monotonic_clock
 *
 * with the license:
 *
 * Copyright (c) 2010 Thomas Habets. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <bfam_base.h>
#include <bfam_log.h>
#include <bfam_clock.h>

#if defined(BFAM_HAVE_CLOCK_MONOTONIC)
int bfam_clock_is_monotonic()
{
  struct timespec ts;

  if (clock_gettime(CLOCK_MONOTONIC, &ts))
  {
    return 0;
  }
  return 1;
}
#elif defined(BFAM_HAVE_MACH_ABSOLUTE_TIME)
int bfam_clock_is_monotonic()
{
  mach_timebase_info_data_t info;

  if (mach_timebase_info(&info))
  {
    return 0;
  }
  return 1;
}
#else
int bfam_clock_is_monotonic() { return 0; }
#endif

static double bfam_clock_fallback()
{
#ifdef BFAM_HAVE_GETTIMEOFDAY
  struct timeval tv;
  static double ref_time_sec = 0;

  if (0 == gettimeofday(&tv, NULL))
  {
    if (ref_time_sec == 0)
    {
      ref_time_sec = (double)tv.tv_sec;
    }
    return (tv.tv_sec - ref_time_sec) + tv.tv_usec / 1000000.0;
  }
  BFAM_WARNING("gettimeofday() failed: %s\n", strerror(errno));
#endif
  return (double)time(0);
}

#if defined(BFAM_HAVE_CLOCK_MONOTONIC)
double bfam_clock()
{
  struct timespec ts;
  static double ref_time_sec = 0;

  if (0 == clock_gettime(CLOCK_MONOTONIC, &ts))
  {
    if (0 == ref_time_sec)
    {
      ref_time_sec = (double)ts.tv_sec;
    }
    return (ts.tv_sec - ref_time_sec) + ts.tv_nsec / 1000000000.0;
  }

  BFAM_WARNING("clock_gettime(CLOCK_MONOTONIC, &ts) failed");
  return bfam_clock_fallback();
}
#elif defined(BFAM_HAVE_MACH_ABSOLUTE_TIME)
double bfam_clock()
{
  uint64_t time = mach_absolute_time();
  static double scaling_factor = 0;
  static uint64_t ref_time = 0;

  if (0 == ref_time)
  {
    ref_time = time;
  }

  time -= ref_time;

  if (scaling_factor == 0)
  {
    mach_timebase_info_data_t info;
    kern_return_t ret = mach_timebase_info(&info);
    if (ret != 0)
    {
      BFAM_WARNING("mach_timebase_info() failed: %d", ret);
      return bfam_clock_fallback();
    }
    scaling_factor = (double)info.numer / (double)info.denom;
  }
  return (double)time * scaling_factor / 1000000000.0;
}
#else
double bfam_clock() { return bfam_clock_fallback(); }
#endif
