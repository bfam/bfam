#ifndef BFAM_CLOCK_H
#define BFAM_CLOCK_H

/** Returns true if the clock time will be monotonic.
 *
 * \return true if clock is monotonic
 */
int bfam_clock_is_monotonic();

/** Returns a clock value that can be used for timing.
 *
 * \return clock value
 */
double bfam_clock();

#endif
