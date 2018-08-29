/*************************************
 * Header for building without CMAKE *
 *************************************/

/* You have to
 * 1. set defines.
 * 2. add header folder to project path
 *    or keep all sources and headers together.
 * 3. link libraries you use manually.
 * */

/**** Defines ****/
/* set 1 for your OS
 * MAC is considered as UNIX
 */
#define OS_UNIX   0
#define OS_WIN    0

/* set 1 if you use BLAS
 * */
#define USE_CBLAS 0

/* set 1 for your BLAS
 * */
#define USE_OPEN  0
#define USE_MKL   0

