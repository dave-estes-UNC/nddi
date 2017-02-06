//
//  FrameVolume.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_FrameVolume_h
#define pixelbridge_FrameVolume_h

#include <cstdlib>
#include <cassert>

#include "NDimensionalDisplayInterface.h"

using namespace std;

namespace nddi {

    class FrameVolume {

    protected:
        CostModel     * costModel_;
        unsigned int    dimensionality_;
        unsigned int  * dimensionalSizes_;
        size_t          size_;
        Pixel         * pixels_;

    public:

        FrameVolume(CostModel* costModel,
                    unsigned int frameVolumeDimensionality,
                    unsigned int* frameVolumeDimensionalSizes)
        : costModel_(costModel), size_(1), pixels_(NULL) {

            dimensionality_ = frameVolumeDimensionality;
            dimensionalSizes_ = (unsigned int*)malloc(sizeof(unsigned int) * frameVolumeDimensionality);
            for (int i = 0; i < dimensionality_; i++) {
                dimensionalSizes_[i] = frameVolumeDimensionalSizes[i];
                size_ *= dimensionalSizes_[i];
            }
            if (!costModel_->isHeadless()) {
                pixels_ = (Pixel *)malloc(sizeof(Pixel) * size_);
                memset(pixels_, 0x00, sizeof(Pixel) * size_);
            }
        }

        ~FrameVolume() {

            if (dimensionalSizes_) {
                free((void*)dimensionalSizes_);
            }
            if (pixels_) {
                free((void*)pixels_);
            }
        }

        unsigned int getSize() {

            return size_;
        }

        void PutPixel(Pixel p, unsigned int* location) {
            if (!costModel_->isHeadless()) {
                setPixel(location, p);
            } else {
                costModel_->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                                     1,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     1 * BYTES_PER_PIXEL,
                                                     0);
            }
        }

        void CopyPixelStrip(Pixel* p, unsigned int* start, unsigned int* end) {
            int dimensionToCopyAlong;
            bool dimensionFound = false;

            // Find the dimension to copy along
            for (int i = 0; !dimensionFound && (i < dimensionality_); i++) {
                if (start[i] != end[i]) {
                    dimensionToCopyAlong = i;
                    dimensionFound = true;
                }
            }

            if (!costModel_->isHeadless()) {
                unsigned int position[dimensionality_];
                memcpy(position, start, sizeof(unsigned int) * dimensionality_);
                for (int j = 0; j <= end[dimensionToCopyAlong] - start[dimensionToCopyAlong]; j++) {
                    setPixel(position, p[j]);
                    position[dimensionToCopyAlong]++;
                }
            } else {
                costModel_->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                                     end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     (end[dimensionToCopyAlong] - start[dimensionToCopyAlong] + 1) * BYTES_PER_PIXEL,
                                                     0);
            }
        }

        void CopyPixels(Pixel* p, unsigned int* start, unsigned int* end) {

            unsigned int position[dimensionality_];
            bool copyFinished = false;
            int pixelsCopied = 0;

            memcpy(position, start, sizeof(unsigned int) * dimensionality_);

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set pixel in frame volume at position
                if (!costModel_->isHeadless())
                    setPixel(position, p[pixelsCopied]);
                pixelsCopied++;

                // Move to the next position
                int fvDim = 0;
                bool overflow;
                do {
                    overflow = false;
                    position[fvDim]++;
                    if ( (position[fvDim] >= dimensionalSizes_[fvDim])
                        || (position[fvDim] > end[fvDim]) ) {
                        overflow = true;
                        position[fvDim] = start[fvDim];
                        if (++fvDim >= dimensionality_)
                            copyFinished = true;
                    }
                } while (overflow && !copyFinished);

            } while (!copyFinished);

            if (costModel_->isHeadless())
                costModel_->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                                     pixelsCopied,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     pixelsCopied * BYTES_PER_PIXEL,
                                                     0);
        }

        void FillPixel(Pixel p, unsigned int* start, unsigned int* end) {

            unsigned int position[dimensionality_];
            bool fillFinished = false;
            int pixelsFilled = 0;

            memcpy(position, start, sizeof(unsigned int) * dimensionality_);

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set pixel in frame volume at position
                if (!costModel_->isHeadless())
                    setPixel(position, p);
                pixelsFilled++;

                // Move to the next position
                int fvDim = 0;
                bool overflow;
                do {
                    overflow = false;
                    position[fvDim]++;
                    if ( (position[fvDim] >= dimensionalSizes_[fvDim])
                        || (position[fvDim] > end[fvDim]) ) {
                        overflow = true;
                        position[fvDim] = start[fvDim];
                        if (++fvDim >= dimensionality_)
                            fillFinished = true;
                    }
                } while (overflow && !fillFinished);

            } while (!fillFinished);

            if (costModel_->isHeadless())
                costModel_->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                                     pixelsFilled,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     pixelsFilled * BYTES_PER_PIXEL,
                                                     0);

        }

        void CopyFrameVolume(unsigned int* start, unsigned int* end, unsigned int* dest) {

            unsigned int positionFrom[dimensionality_];
            unsigned int positionTo[dimensionality_];
            bool copyFinished = false;
            int pixelsCopied = 0;

            memcpy(positionFrom, start, sizeof(unsigned int) * dimensionality_);
            memcpy(positionTo, dest, sizeof(unsigned int) * dimensionality_);

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set pixel in frame volume at position
                if (!costModel_->isHeadless())
                    setPixel(positionTo, getPixel(positionFrom));
                pixelsCopied++;

                // Move to the next position
                int fvDim = 0;
                bool overflow;
                do {
                    overflow = false;
                    positionFrom[fvDim]++;
                    positionTo[fvDim]++;
                    if ( (positionFrom[fvDim] >= dimensionalSizes_[fvDim])
                        || (positionFrom[fvDim] > end[fvDim]) ) {
                        overflow = true;
                        positionFrom[fvDim] = start[fvDim];
                        positionTo[fvDim] = dest[fvDim];
                        if (++fvDim >= dimensionality_)
                            copyFinished = true;
                    }
                } while (overflow && !copyFinished);

            } while (!copyFinished);

            if (costModel_->isHeadless())
                costModel_->registerBulkMemoryCharge(FRAME_VOLUME_COMPONENT,
                                                     pixelsCopied,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     pixelsCopied * BYTES_PER_PIXEL,
                                                     0);

        }

        void setPixel(unsigned int* location, Pixel pixel) {

            unsigned int offset = 0;
            unsigned int multiplier = 1;

            assert(!costModel_->isHeadless());

            for (int i = 0; i < dimensionality_; i++) {
                assert(location[i] < dimensionalSizes_[i]);

                offset += location[i] * multiplier;
                multiplier *= dimensionalSizes_[i];
            }

            pixels_[offset].packed = pixel.packed;

            costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, WRITE_ACCESS, pixels_ + offset, BYTES_PER_PIXEL, 0);
        }

        Pixel getPixel(unsigned int* location) {

            Pixel         pixel;
            unsigned int  offset = 0;
            unsigned int  multiplier = 1;

            assert(!costModel_->isHeadless());

            for (int i = 0; i < dimensionality_; i++) {
                assert(location[i] < dimensionalSizes_[i]);

                offset += location[i] * multiplier;
                multiplier *= dimensionalSizes_[i];
            }

            pixel.packed = pixels_[offset].packed;

            costModel_->registerMemoryCharge(FRAME_VOLUME_COMPONENT, READ_ACCESS, pixels_ + offset, BYTES_PER_PIXEL, 0);

            return pixel;
        }

        Pixel * data() {

            return pixels_;
        }
    };
}

#endif
