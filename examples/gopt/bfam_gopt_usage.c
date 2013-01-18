#include <bfam.h>

int
main(int argc, const char **argv)
{

  FILE *out_stream;
  const char *filename;
  int verbosity;
  int i;
  void *options= bfam_gopt_sort( & argc, argv, bfam_gopt_start(
      bfam_gopt_option( 'h', 0, bfam_gopt_shorts( 'h', '?' ), bfam_gopt_longs( "help", "HELP" )),
      bfam_gopt_option( 'z', 0, bfam_gopt_shorts( 'V' ), bfam_gopt_longs( "version" )),
      bfam_gopt_option( 'v', BFAM_GOPT_REPEAT, bfam_gopt_shorts( 'v' ), bfam_gopt_longs( "verbose" )),
      bfam_gopt_option( 'o', BFAM_GOPT_ARG, bfam_gopt_shorts( 'o' ), bfam_gopt_longs( "output" ))));
  /*
   * there are four possible options to this program, some of which have
   * multiple names:
   *
   * -h -? --help --HELP
   * -V --version
   * -v --verbose  (which may be repeated for more verbosity)
   * -o --output  (which requires an option argument)
   *
   * the program will have been terminated if unrecognised options are
   * specified, or if an option without BFAM_GOPT_REPEAT was repeated, or if an
   * option without BFAM_GOPT_ARG had an option argument, or if one with
   * BFAM_GOPT_ARG had no argument.
   */

  if( bfam_gopt( options, 'h' ) ){
    /*
     * if any of the help options was specified
     */
    fprintf( stdout, "help text\n" );
    exit( EXIT_SUCCESS );
  }

  if( bfam_gopt( options, 'z' ) ){
    /*
     * if --version was specified
     * NB: 'z' is just a key within the source code
     * if you specify -z, it will be treated as an unknown option
     */
    fprintf( stdout, "version number\n" );
    exit( EXIT_SUCCESS );
  }

  if( bfam_gopt_arg( options, 'o', & filename ) && strcmp( filename, "-" ) ){
    /*
     * if -o or --output was specified, and its argument was not "-"
     */
    out_stream= fopen( filename, "wb" );
    if( ! out_stream ){
      fprintf( stderr, "%s: %s: could not open file for output\n", argv[0],
          filename );
      exit( EXIT_FAILURE);
    }
  }
  else
    out_stream= stdout;

  verbosity= bfam_gopt( options, 'v' );
  /*
   * return value is the number of times that the option was specified
   */

  if( verbosity > 1 )
    fprintf( stderr, "being really verbose\n" );

  else if( verbosity )
    fprintf( stderr, "being verbose\n" );

  bfam_gopt_free( options );
  /*
   * release memory used
   * you can no longer call bfam_gopt() etc.
   */

  for( i= 0; i < argc; ++i )
    fprintf( out_stream, "%s\n", argv[i] );
    /*
     * all the options have been removed
     * leaving only the program name and operands
     */

  exit( EXIT_SUCCESS );
}
/* To accept an option that can be repeated and takes an argument, use
 * BFAM_GOPT_REPEAT | BFAM_GOPT_ARG.  If you do this then you will want to use
 * bfam_gopt_arg_i() or bfam_gopt_args().  See bfam_gopt.h for details.
 */
