#ifndef NDDIFEATURES_H_
#define NDDIFEATURES_H_

/*
 * These features are all on by default. Uncomment a line (and cross fingers) to turn one off.
 */
#ifdef DEBUG
#define NO_OMP
#endif

//#define NO_OMP
//#define NO_GL
#define NO_CL


/*
 * Running make with NO_HACKS=1 should turn these off if using the Makefile for Linux,
 * but uncommenting these will accomplish the same thing. And in OS X it's the only way
 * since these hacks are hardcoded in the project for lack of recent Makefile support
 * for OS X.
 */
//#undef SUPRESS_EXCESS_RENDERING
//#undef SKIP_COMPUTE_WHEN_SCALER_ZERO

/*
 * This turns on the OpenCL profiling code which updates the cost model with time
 * NOTE: This triggers a lot of clFinish()
 */
//#define CL_PROFILING_ENABLED

/*
 * Turns on support for alpha channel support for coefficient plane blending
 */
#define USE_ALPHA_CHANNEL

/*
 * Use the narrowed data fields instead of 4-byte words for EVERYTHING.
 */
#define USE_NARROW_DATA_FIELDS

/*
 * Uses the nDDI extension to update groups of tiles
 */
#define USE_COPY_PIXEL_TILES

/*
 * Used to dramatically narrow the various data stores. Can lead to bugs, so proceed carefully.
 */
#define NARROW_DATA_STORES

#endif /* NDDIFEATURES_H_ */