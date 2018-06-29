#ifndef NDDIFEATURES_H_
#define NDDIFEATURES_H_

// Debugging is very hard with OpenMP enabled.
#ifdef DEBUG
#undef USE_OMP
#endif

//#define NO_GL
#define NO_CL


/*
 * Running cmake with -DHACKS=off ... should turn these off, but uncommenting
 * these will accomplish the same thing individually.
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
 * For cost model calculations, use the narrowed data fields instead of 4-byte words for
 * EVERYTHING. This doesn't affect any data stores or calculations.
 */
#define COST_MODEL_USE_NARROW_DATA_FIELDS

/*
 * Uses the nDDI extension to update groups of tiles
 */
#define USE_COPY_PIXEL_TILES

/*
 * Used to dramatically narrow the various data stores. Can lead to bugs, so proceed
 * carefully. For instance, coefficients are 4 bytes by design, but this flag will
 * reduce them to two bytes. If the use case needs a coefficient beyond 2 bytes, then
 * there's likely going to be an unguarded overflow.
 */
#define NARROW_DATA_STORES

#endif /* NDDIFEATURES_H_ */
