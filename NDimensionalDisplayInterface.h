#ifndef N_DIMENSIONAL_DISPLAY_INTERFACE_H
#define N_DIMENSIONAL_DISPLAY_INTERFACE_H

/**
 * \file NDimensionalDisplayInterface.h
 *
 * \brief This file embodies the nDDI C++ API.
 *
 * This file embodies the nDDI C++ API. Any C++ implementation must extend this class and use the
 * data types, enumerations, and macros within. Implementations might range from the user space API
 * with an nDDI kernel device driver to a networked nDDI client communicating with a remote nDDI
 * server.
 */

#include <vector>
#include <stdint.h>
#include "CostModel.h"
///\cond
#include <cereal/archives/xml.hpp>
///\endcond

using namespace std;

/**
 * \brief Namespace for the entire nDDI API.
 *
 * Namespace for the entire nDDI API.
 */
namespace nddi {

    /**
     * \brief Special value when passing coefficients indicating that the coefficient in that location should
     *        remain unchanged.
     *
     * When specifying coefficient matrices for purposes of updating the coefficient plane, the NDDI client
     * can use this value for one or more of the elements in the matrix if they would like the element in the
     * same location of the destination coefficient matrix to remain unchanged.
     */
    #define COEFFICIENT_UNCHANGED INT16_MAX

    /**
     * \brief Special value when passing coefficients indicating that the coefficient's X coordinate should
     *        be used as its value.
     *
     * Special value when passing coefficients indicating that the coefficient's X coordinate should
     * be used as its value.
     */
    #define COEFFICIENT_MATRIX_X (INT16_MIN + 2)

    /**
     * \brief Special value when passing coefficients indicating that the coefficient's Y coordinate should
     *        be used as its value.
     *
     * Special value when passing coefficients indicating that the coefficient's Y coordinate should
     * be used as its value.
     */
    #define COEFFICIENT_MATRIX_Y (INT16_MIN + 1)
  
    /**
     * \brief Special value when passing coefficients indicating that the coefficient's P (plane) coordinate should
     *        be used as its value.
     *
     * Special value when passing coefficients indicating that the coefficient's P (plane) coordinate should
     * be used as its value.
     */
    #define COEFFICIENT_MATRIX_P (INT16_MIN + 0)

    /**
     * \brief The default value indicating a full scaler.
     *
     * The default value indicated a full scaler. Can be changed by the client when configuring the
     * nDDI display.
     */
    #define DEFAULT_FULL_SCALER 256

    /**
     * \brief Union representing an RGBA 32-bit pixel.
     *
     * Union representing an RGBA 32-bit pixel.
     * Implementation may drop the alpha channel.
     */
    typedef union {
        /**
        * \brief Struct used to access individual bytes of the 32-bit pixel.
        *
        * Struct used to access individual bytes of the 32-bit pixel.
        */
        struct {
            /// \brief Red field
            ///
            /// Red field
            uint8_t r;
            /// \brief Green field
            ///
            /// Green field
            uint8_t g;
            /// \brief Blue field
            ///
            /// Blue field
            uint8_t b;
            /// \brief Alpha field
            ///
            /// Alpha field
            uint8_t a;
        } fields;
      
        /**
         * \brief 32-bit wide field to access the entire pixel as a 32-bit value.
         *
         * 32-bit wide field to access the entire pixel as a 32-bit value.
         */
        uint32_t packed;

        ///\cond
        template <class Archive>
        void serialize(Archive& ar) {
          ar(CEREAL_NVP(r), CEREAL_NVP(g), CEREAL_NVP(b), CEREAL_NVP(a));
        }
        ///\endcond
    } Pixel;

    /**
     * \brief Union representing the 4-channel scaler.
     *
     * Union representing the 4-channel scaler.
     * Implementation may drop the alpha channel.
     */
    typedef union {
        /**
         * \brief Struct used to access individual 16-bit shorts of the 64-bit scaler.
         *
         * Struct used to access individual 16-bit shorts of the 64-bit scaler.
         */
        struct {
            /// \brief Red field
            ///
            /// Red field
            int16_t r;
            /// \brief Green field
            ///
            /// Green field
            int16_t g;
            /// \brief Blue field
            ///
            /// Blue field
            int16_t b;
            /// \brief Alpha field
            ///
            /// Alpha field
            int16_t a;
        } fields;
      
        /**
        * \brief 64-bit wide field to access the entire pixel as a 64-bit value.
        *
        * 64-bit wide field to access the entire pixel as a 64-bit value.
        */
        uint64_t packed;

        ///\cond
        template <class Archive>
        void serialize(Archive& ar) {
          ar(CEREAL_NVP(r), CEREAL_NVP(g), CEREAL_NVP(b), CEREAL_NVP(a));
        }
        ///\endcond
    } Scaler;

    /**
     * \brief Options for the sign mode of the pixel byte channels.
     *
     * Options for the sign mode of the pixel byte channels. Configuring an nDDI display to use signed
     * pixels treats each pixel channel as a signed 8-bit value when multiplying and
     * accumulating a pixel value across all of the Coefficient Planes. Default is
     * UNSIGNED_MODE.
     */
    typedef enum {
        UNSIGNED_MODE = 0,
        SIGNED_MODE = 1
    } SignMode;

    /**
     * \brief Type for the coefficients stored in the coefficient matrices.
     *
     * Type for the coefficients stored in the coefficient matrices.
     */
#ifdef NARROW_DATA_STORES
    typedef int16_t Coeff;
#else
    typedef int32_t Coeff;
#endif

      
    /**
     * \brief This abstract class serves as a software interface to an n-Dimensional Display Interface (nDDI) compliant
     *        display device.
     *
     * This abstract class serves as a software interface to an n-Dimensional Display Interface (nDDI) compliant
     * display device. Implementations of this interface may work with a simulated display that writes
     * to a system framebuffer, an embedded display that works with a device driver, or a remote display
     * connected via an IP-based socket and communicating with an NDDI specific protocol.
     *
     */
    class NDimensionalDisplayInterface  {

    public:

        /**
         * \brief Required default constructor for abstract class NDimensionalDisplayInterface.
         *
         * Required default constructor for abstract class NDimensionalDisplayInterface.
         */
        NDimensionalDisplayInterface() {}

        /**
         * \brief Minimal constructor which uses a fixed display size.
         *
         * Minimal constructor which uses a fixed display size and configures the nDDI display with the provided
         * dimensions for the Input Vector, Coefficient Planes, and Frame Volume.
         *
         * @param frameVolumeDimensionalSizes This vector is used to configure the frame volume.
         *                                    Each element in the vector represents a dimension and that element's
         *                                    value represents the size of that dimension. e.g. a simple 4x4 2D
         *                                    frame volume will be configured with a two-element vector with 4 and 4 in it.
         * @param numCoefficientPlanes Specifies how many of the maximum coefficient planes will be used.
         * @param inputVectorSize Used to configure the size of the input vector. It must be greater than or equal to two.
         * @param fixed8x8MacroBlocks Configuration option which can configure the nDDI display to have fixed macroblocks.
         *                            This is strictly used for memory optimizations so that one coefficient matrix is used
         *                            for an entire macro block instead of one per pixel location. Note: This is still per plane.
         * @param useSingleCoefficientPlane Configuration option which can configure the display implementation to use only one
         *                                  Coefficient Plane even when it's configured to used multiple planes. The configured
         *                                  number of planes are still all logically used for blending, but their coefficients may
         *                                  be different if using special values COEFFICIENT_MATRIX_X, COEFFICIENT_MATRIX_Y, or
         *                                  COEFFICIENT_MATRIX_P. This is strictly use for memory optimizations in simulation.
         */
        NDimensionalDisplayInterface(vector<unsigned int> &frameVolumeDimensionalSizes,
                                     unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                                     bool fixed8x8Macroblocks, bool useSingleCoeffcientPlane) {}
        /**
         * \brief Full constructor which additionally allows the display size to be configured.
         *
         * Full constructor which additionally allows the display size to be configured.
         *
         * @param frameVolumeDimensionalSizes This vector is used to configure the frame volume.
         *                                    Each element in the vector represents a dimension and that element's
         *                                    value represents the size of that dimension. e.g. a simple 4x4 2D
         *                                    frame volume will be configured with a two-element vector with 4 and 4 in it.
         * @param displayWidth Used to configure the width of the display if it is less than the display device.
         * @param displayHeight Used to configure the width of the display if it is less than the display device.
         * @param numCoefficientPlanes Specifies how many of the maximum coefficient planes will be used.
         * @param inputVectorSize Used to configure the size of the input vector. It must be greater than or equal to two.
         * @param fixed8x8MacroBlocks Configuration option which can configure the nDDI display to have fixed macroblocks.
         *                            This is strictly used for memory optimizations so that one coefficient matrix is used
         *                            for an entire macro block instead of one per pixel location. Note: This is still per plane.
         * @param useSingleCoefficientPlane Configuration option which can configure the display implementation to use only one
         *                                  Coefficient Plane even when it's configured to used multiple planes. The configured
         *                                  number of planes are still all logically used for blending, but their coefficients may
         *                                  be different if using special values COEFFICIENT_MATRIX_X, COEFFICIENT_MATRIX_Y, or
         *                                  COEFFICIENT_MATRIX_P. This is strictly use for memory optimizations in simulation.
         */
        NDimensionalDisplayInterface(vector<unsigned int> &frameVolumeDimensionalSizes,
                                     unsigned int displayWidth, unsigned int displayHeight,
                                     unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                                     bool fixed8x8Macroblocks, bool useSingleCoeffcientPlane) {}

        /**
         * \brief Used to query the display width.
         *
         * Used to query the display width.
         *
         * @return The width of the display.
         */
        virtual unsigned int DisplayWidth() = 0;

        /**
         * \brief Used to query the display height.
         *
         * Used to query the display height.
         *
         * @return The height of the display.
         */
        virtual unsigned int DisplayHeight() = 0;

        /**
         * \brief Used to query the number of coefficient planes.
         *
         * Used to query the number of coefficient planes.
         *
         * @return The number of coefficient planes.
         */
        virtual unsigned int NumCoefficientPlanes() = 0;

        /**
         * \brief Copies the provided pixel to the specified location.
         *
         * Copies the provided pixel to the specified location.
         *
         * @param p The pixel value to be copied.
         * @param location Tuple for the location within the frame volume where the pixel will be copied to.
         */
        virtual void PutPixel(Pixel p, vector<unsigned int> &location) = 0;

        /**
         * \brief Copies the one dimensional array of pixels along a particular dimension in the frame volume.
         *
         * Copies the one dimensional array of pixels along a particular dimension in the frame volume. In a
         * two-dimensional frame volume, this can be thought of as a way to copy along a row or along a column, but
         * not both since the input pixels are only one-dimensional.
         *
         * @param p The pointer to the pixel values to be copied.
         * @param start Tuple for the first pixel in the frame volume to be filled.
         * @param end Tuple for the last pixel in the frame volume to be filled. All but one of the values in
         *            values in this last pixel should be identical to the start pixel.
         */
        virtual void CopyPixelStrip(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * \brief Copies the array of pixels into the designated region of the frame volume.
         *
         * Copies the array of pixels into the designated region of the frame volume. The data must be
         * arranged in the array with strides for each dimension of the area. So to copy pixels into a
         * 2 x 2 x 2 region in the frame volume, the array must be arranged accordingly:
         * (0,0,0) (1,0,0) (0,1,0) (1,1,0) (0,0,1) (1,0,1) (0,1,1) (1,1,1)
         *
         * @param p The pointer to the pixel values to be copied.
         * @param start Tuple for the first pixel in the frame volume to be filled.
         * @param end Tuple for the last pixel in the frame volume to be filled.
         */
        virtual void CopyPixels(Pixel* p, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * \brief Copies the array of pixels into the designated tile regions of the frame volume.
         *
         * Copies the array of pixels into the designated tile regions of the frame volume. The data must be
         * arranged in the array with strides for each dimension of the area. Only 2D tiles are supported.
         *
         * @param p The pointer to the pixel values to be copied.
         * @param starts Vector holding series of tuples for the first pixel for each destination tile in the frame volume.
         * @param size Two element tuple for the size of each tile (w, h).
         */
        virtual void CopyPixelTiles(vector<Pixel*> &p, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) = 0;

        /**
         * \brief Fills the frame volume with the specified pixel.
         *
         * Fills the frame volume with the specified pixel. It can fill in multiple
         * dimensions by starting at the start pixel and filling in each dimension until
         * the end pixel value is reached.
         *
         * @param p The pixel value to be filled.
         * @param start Tuple for the first pixel in the frame volume to be filled.
         * @param end Tuple for the last pixel in the frame volume to be filled.
         */
        virtual void FillPixel(Pixel p, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * \brief Copies pixels from one multi-dimensional region of the frame volume to another region.
         *
         * Copies pixels from one multi-dimensional region of the frame volume to another region.
         *
         * @param start Tuple for the starting coordinate of the source region.
         * @param end Tuple for the ending coordinate of the source region.
         * @param dest Tuple for the first starting pixel of the destination region to be filled.
         */
        virtual void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest) = 0;

        /**
         * \brief Used to update the input vector with the extra values in the input vector.
         *
         * Used to update the input vector with the extra values in the input vector.
         *
         * @param input Tuple for the values to use for the update. The length of this
         *              tuple must equal the size of the actual input vector
         *              minus two, since the first two values in the input
         *              vector cannot be changed.
         */
        virtual void UpdateInputVector(vector<int> &input) = 0;

        /**
         * \brief Used to copy the specified coefficientMatrix into the specified location of the coefficient
         *        planes.
         *
         * Used to copy the specified coefficientMatrix into the specified location of the coefficient
         * planes.
         *
         * @param coefficientMatrix This two-dimensional vector holds the matrix to be copied.
         *                          It's size must match the configuration of the coefficient matrices
         *                          exactly. Can use COFFICIENT_UNCHANGED for one or more elements.
         * @param location This two-element vector specifies the tuple for the location in the coefficient plane where the provided
         *                 coefficient matrix will be copied.
         */
        virtual void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) = 0;

        /**
         * \brief Used to copy the specified coefficientMatrix into a range of locations in the coefficient planes.
         *
         * Used to copy the specified coefficientMatrix into a range of locations in the coefficient planes.
         *
         * @param coefficientMatrix This two-dimensional vector holds the matrix to be copied.
         *                          It's size must match the configuration of the coefficient matrices
         *                          exactly. Can use COFFICIENT_UNCHANGED for one or more elements.
         * @param start This three-element vector specifies the tuple of location in the coefficient planes where the first
         *              coefficient matrix will be copied to.
         * @param end This three-element vector specifies the tuple of the location in the coefficient planes where the last
         *            coefficient matrix will be copied to.
         */
        virtual void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * \brief Used to copy the specified single coefficient value from a matrix into a range of locations in the coefficient planes.
         *
         * Used to copy the specified single coefficient value from a matrix into a range of locations in the coefficient planes.
         *
         * @param coefficient This single value will be placed into each coefficient at the specified location in the coefficient
         *                    matrices of the specified range.
         * @param row The row of the coefficient to be updated in the coefficient matrix.
         * @param col The column of the coefficient to be updated in the coefficient matrix.
         * @param start This three-element vector specifies the tuple of location in the coefficient planes where the first
         *              coefficient matrix will be copied to.
         * @param end This three-element vector specifies the tuple of location in the coefficient planes where the last
         *            coefficient matrix will be copied to.
         */
        virtual void FillCoefficient(int coefficient, unsigned int row, unsigned int col, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * \brief For each coefficient, positions, and start; copies the coefficient to the position
         *        in the in each coefficient matrix in the 2D tile specified by the start and size.
         *
         * For each coefficient, positions, and start; copies the coefficient to the position
         * in the in each coefficient matrix in the 2D tile specified by the start and size.
         *
         * @param coefficients The buffer of coefficients.
         * @param positions Tuple for the position (row, col) to place the coefficient within the coefficient matrix.
         * @param starts Tuple for the location (x, y) of the start of the tile in the coefficient planes.
         * @param size Tuple for the size (w, h) of the tile.
         */
        virtual void FillCoefficientTiles(vector<int> &coefficients, vector<vector<unsigned int> > &positions, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) = 0;

        /**
         * \brief Used to copy the specified scaler to a range of locations in the coefficient planes.
         *
         * Used to copy the specified scaler to a range of locations in the coefficient planes.
         *
         * @param scaler This single scaler will be copied to each location in the range of coefficient planes.
         * @param start This three-element vector specifies tuple for the start of the range in the coefficient planes
         *              where the scaler will be copied to.
         * @param end This three-element vector specifies the tuple for end of the range in the coefficient planes where
         *            the scalers will be copied to.
         */
        virtual void FillScaler(Scaler scaler, vector<unsigned int> &start, vector<unsigned int> &end) = 0;

        /**
         * \brief Used to copy the specified scalers to a series of 2D ranges of locations (tiles) in the coefficient planes.
         *
         * Used to copy the specified scalers to a series of 2D ranges of locations (tiles) in the coefficient planes.
         * This is accomplished with a set of scalers, an equal number of tile starts, and one tile size.
         *
         * @param scalers Each scaler in this list will be filled to its own tile, which is a 2D range of coefficient matrices.
         * @param starts Vector of tuples for the start locations (x, y, z) in the coefficient planes for each tile to be filled.
         * @param size Tuple for the size (w, h) of the tile.
         */
        virtual void FillScalerTiles(vector<uint64_t> &scalers, vector<vector<unsigned int> > &starts, vector<unsigned int> &size) = 0;

        /**
         * \brief Used to copy the specified scalers to a stack of 2D ranges of locations (tiles) in the coefficient planes.
         *
         * Used to copy the specified scalers to a stack of 2D ranges of locations (tiles) in the coefficient planes.
         * This is accomplished with with a set of scalers, a single tile stack location, and one tile size.
         * The stack includes the top-most tile on the coefficient plane indicated by the tile start as well as
         * the tiles for the planes under (higher plane coordinates) the start tile. The height of the stack to be
         * filled is determined by the number of scalers provided.
         *
         * @param scalers Each scaler in this list will be filled to its own tile (2D range of coefficient matrices) in the stack.
         * @param start Tuple for the location (x, y) of the start of the tile stack in the coefficient planes.
         * @param size Tuple for the size (w, h) of the tile.
         */
        virtual void FillScalerTileStack(vector<uint64_t> &scalers, vector<unsigned int> &start, vector<unsigned int> &size) = 0;

        /**
         * \brief Allows the bytes of pixel values to be interpretted as signed values when scaling, accumulating, and clamping
         * in the pixel blending pipeline.
         *
         * Allows the bytes of pixel values to be interpretted as signed values when scaling, accumulating, and clamping
         * in the pixel blending pipeline. When the SignMode is set to the default UNSIGNED_MODE, each 8-bit color channel for
         * pixels will be treated as an unsigned value. Setting this configuration option to SIGNED_MODE treats those channels
         * as 8-bit signed values.
         *
         * @param mode Can be UNSIGNED_MODE or SIGNED_MODE.
         */
        virtual void SetPixelByteSignMode(SignMode mode) = 0;

        /**
         * \brief Used to set the scaler value which is interpretted as fully on or 100%.
         *
         * Used to set the scaler value which is interpretted as fully on or 100%. The default is
         * 256, which implies that any scaler sent by the client is an
         * integer fraction of 256, but in fact a scaler can be larger than 256,
         * leading to planes that contribute 2.5x or even -3x for instance.
         *
         * @param scaler The value to be interpretted as fully on or 100%.
         */
        virtual void SetFullScaler(uint16_t scaler) = 0;

        /**
         * \brief Used to get the current full scaler value.
         *
         * Used to get the current full scaler value.
         *
         * @return The current fully on scaler value.
         */
        virtual uint16_t GetFullScaler() = 0;

        /**
         * \brief Returns the CostModel for this display.
         *
         * Returns the CostModel for this display. The CostModel can be queried by the
         * host application to understand the cost of operations after they complete.
         * The CostModel is purely a mechanism for collecting experiment data for
         * simulated nDDI displays.
         *
         * @return The CostModel created and maintained by this display.
         */
        virtual CostModel* GetCostModel() = 0;
    };

} // namespace nddi {

#endif // N_DIMENSIONAL_DISPLAY_INTERFACE_H
