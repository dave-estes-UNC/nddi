#ifndef N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H
#define N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H

/**
 * \file NDimensionalDisplayInterfaceExtended.h
 *
 * \brief This file embodies the extensions to the nDDI C++ API.
 *
 * This file embodies the extensions to the nDDI C++ API.
 */

/**
 * \brief Namespace for the entire nDDI API.
 *
 * Namespace for the entire nDDI API.
 */
namespace nddi {
    
    /**
     * \brief This abstract class serves as an extension to the NDDI software interface.
     *
     * This abstract class serves as an extension to the NDDI software interface.
     */
    class NDimensionalDisplayInterfaceExtended  {
        
    public:

        /**
         * \brief Copies from one region of the frame volume to another with blending.
	 *
	 * Copies from one region of the frame volume to another with blending.
         *
         * @param start The first pixel in the frame volume to be copied from.
         * @param end The last pixel in the frame volume to be copied from.
         * @param dest The first pixel in the frame volume to be copied to.
         * @param blend If true, then the copy then the alpha channel for each src
	 *              pixel will be used to blend the source range over the
	 *              destination range.
         */
    	virtual void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest, bool blend) = 0;
};
    
} // namespace nddi {

#endif // N_DIMENSIONAL_DISPLAY_INTERFACE_EXTENDED_H
