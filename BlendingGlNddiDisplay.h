//
//  BlendingGlNddiDisplay.h
//  pixelbridge
//
//  Created by Dave Estes on 1/31/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_BlendingGlNddiDisplay_h
#define pixelbridge_BlendingGlNddiDisplay_h

#include "GlNddiDisplay.h"
#include "NDimensionalDisplayInterfaceExtended.h"

using namespace nddi;
using namespace std;

/**
 * Blending version of a GL NDDI Display.
 */
class BlendingGlNddiDisplay : public GlNddiDisplay, public NDimensionalDisplayInterfaceExtended {

public:
    BlendingGlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                          unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                          bool headless = false);
    BlendingGlNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                          unsigned int displayWidth, unsigned int displayHeight,
                          unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                          bool headless = false);
    ~BlendingGlNddiDisplay();

    // To satisfy the NDimensionalDisplayInterfaceExtended interface
    void CopyFrameVolume(vector<unsigned int> &start, vector<unsigned int> &end, vector<unsigned int> &dest, bool blend);

#ifndef SUPRESS_EXCESS_RENDERING
protected:

    /**
     * Renders each pixel of the frame buffer by setting the x, y in the input vector and computing which
     * frame volume pixel to use.
     */
    virtual void Render() = 0;
#endif

private:
    nddi::Pixel BlendPixel(nddi::Pixel pTo, nddi::Pixel pFrom);

};

#endif
