#include <iostream>
#include <stdio.h>
#include <sys/time.h>

#include "Features.h"
#include "SimpleNddiDisplay.h"

inline uint8_t CLAMP_SIGNED_BYTE(int32_t i) {
    uint32_t ret;

    if (i < 0) {
        ret = 0;
    } else if (i > 0xff) {
        ret = 0xff;
    } else {
        ret = i;
    }

    return (uint8_t)ret;
}

inline uint8_t CLAMP_UNSIGNED_BYTE(uint32_t i) {
    uint32_t ret;

    if (i > 0xff) {
        ret = 0xff;
    } else {
        ret = i;
    }

    return ret;
}

inline uint8_t TRUNCATE_BYTE(int32_t i) {
    uint32_t ret;

    ret = i & 0xff;

    return (uint8_t)ret;
}

// public

SimpleNddiDisplay::SimpleNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                             bool headless, unsigned char logcosts, bool fixed8x8Macroblocks, bool useSingleCoeffcientPlane) {
    SimpleNddiDisplay(frameVolumeDimensionalSizes, 320, 240, numCoefficientPlanes, inputVectorSize);
}

SimpleNddiDisplay::SimpleNddiDisplay(vector<unsigned int> &frameVolumeDimensionalSizes,
                             unsigned int displayWidth, unsigned int displayHeight,
                             unsigned int numCoefficientPlanes, unsigned int inputVectorSize,
                             bool headless, unsigned char logcosts, bool fixed8x8Macroblocks, bool useSingleCoeffcientPlane)
: BaseNddiDisplay(frameVolumeDimensionalSizes, displayWidth, displayHeight, numCoefficientPlanes, inputVectorSize, headless, logcosts, fixed8x8Macroblocks, useSingleCoeffcientPlane) {

    numPlanes_ = numCoefficientPlanes;
    frameVolumeDimensionalSizes_ = frameVolumeDimensionalSizes;
    displayWidth_ = displayWidth;
    displayHeight_ = displayHeight;
    pixelSignMode_ = UNSIGNED_MODE;
    quiet_ = true;

    // Create the CostModel
    costModel = new CostModel(frameVolumeDimensionalSizes, inputVectorSize, headless, logcosts);

    // Setup Input Vector
    inputVector_ = new InputVector(costModel, inputVectorSize);

    // Setup framevolume and initialize to black
    frameVolume_ = new FrameVolume(costModel, frameVolumeDimensionalSizes);

    // Setup coefficient plane with zeroed coefficient matrices
    coefficientPlanes_ = new CoefficientPlanes(costModel,
            displayWidth_, displayHeight_, numCoefficientPlanes,
            CM_WIDTH, CM_HEIGHT,
            fixed8x8Macroblocks, useSingleCoeffcientPlane);

    // Setup framebuffer and initialize to black
    frameBuffer_ = (Pixel*)malloc(sizeof(Pixel) * displayWidth_ * displayHeight_);
    memset(frameBuffer_, 0x00, sizeof(Pixel) * displayWidth_ * displayHeight_);

    // Set the full scaler and the accumulator
    SetFullScaler(DEFAULT_FULL_SCALER);

    emptyCoefficientMatrix.resize(inputVectorSize, std::vector<int>(frameVolumeDimensionalSizes.size(), 0));
}

// TODO(CDE): Why is the destructor for SimpleNddiDisplay being called when we're using a ClNddiDisplay?
SimpleNddiDisplay::~SimpleNddiDisplay() {

    delete(inputVector_);
    delete(frameVolume_);
    delete(coefficientPlanes_);

    if (frameBuffer_)
        free(frameBuffer_);
}

// Private

void SimpleNddiDisplay::Render() {
    Render(0, 0, displayWidth_, displayHeight_);
}

void SimpleNddiDisplay::Render(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h) {

    timeval startTime, endTime; // Used for timing data
    if (!quiet_)
        gettimeofday(&startTime, NULL);

    bool doCostCalculation = true;
#ifdef USE_OMP
    doCostCalculation = false;
#pragma omp parallel for
#endif
    for (unsigned int y = sub_y; y < (sub_y + sub_h); y += 8) {
        for (unsigned int x = sub_x; x < (sub_x + sub_w); x += 8) {
            ComputePixels(x, y, 8, doCostCalculation);
        }
    }

    if (!quiet_) {
        gettimeofday(&endTime, NULL);
        printf("Render Statistics:\n  Size: %zdx%zd - FPS: %f\n",
                displayWidth_,
                displayHeight_,
                1.0f / ((double)(endTime.tv_sec * 1000000
                                + endTime.tv_usec
                                - startTime.tv_sec * 1000000
                                - startTime.tv_usec) / 1000000.0f)
                  );
    }
#ifdef SUPRESS_EXCESS_RENDERING
    changed_ = false;
#endif

}

void SimpleNddiDisplay::ComputePixels(unsigned int x, unsigned int y, unsigned int length, bool doCostCalculation) {

    int32_t       rAccumulator, gAccumulator, bAccumulator;
    Pixel         q;
    unsigned int  startx = x, starty = y;
    bool          fixed = fixed8x8Macroblocks_ && length == 8;


    // Current location in the coefficient planes (x, y, p)
    vector<unsigned int> location(3, 0);

    for (y = starty; y < starty + length && y < displayHeight_; y++) {
        location[1] = y;
        for (x = startx; x < startx + length && x < displayWidth_; x++) {
            location[0] = x;
            rAccumulator = gAccumulator = bAccumulator = 0;
            q.packed = 0;

            // Accumulate color channels for the pixels chosen by each plane
            for (unsigned int p = 0; p < numPlanes_; p++) {
                location[2] = p;

                // Grab the scaler for this location
                Scaler scaler;
                scaler = coefficientPlanes_->getScaler(x, y, p);
                if (doCostCalculation) {
                    costModel->registerScalerMemoryCharge(READ_ACCESS, location, location);
                }

#ifdef SKIP_COMPUTE_WHEN_SCALER_ZERO
                if (scaler.packed == 0) continue;
#endif

                // Compute the position vector for the proper pixel in the frame volume.
                vector<unsigned int> fvPosition;
                // Matrix multiply the input vector by the coefficient matrix
                for (int j = 0; j < CM_HEIGHT; j++) {
                    // Initialize to zero
                    fvPosition.push_back(0);
                    // No need to read the x and y from the input vector, just multiply directly.
                    fvPosition[j] += coefficientPlanes_->getCoefficient(location, 0, j) * x;
                    fvPosition[j] += coefficientPlanes_->getCoefficient(location, 1, j) * y;
                    // Then multiply the remainder of the input vector
                    for (int i = 2; i < CM_WIDTH; i++) {
                        fvPosition[j] += coefficientPlanes_->getCoefficient(location, i, j) * inputVector_->getValue(i);
                    }
                }
                if (doCostCalculation) {
                    costModel->registerCoefficientMatrixMemoryCharge(READ_ACCESS, location, location, emptyCoefficientMatrix);
                    if (CM_WIDTH > 2) {
                        for (int j = 0; j < CM_HEIGHT; j++) {
                            costModel->registerInputVectorMemoryCharge(READ_ACCESS, 2, CM_WIDTH - 1);
                        }
                    }
                }

                // Grab the pixel from the frame volume
                q = frameVolume_->getPixel(fvPosition);
                if (doCostCalculation) {
                    costModel->registerFrameVolumeMemoryCharge(READ_ACCESS, fvPosition, fvPosition);
                }

                // Compute the pixel's contribution and add to the accumulator.
        #ifdef USE_ALPHA_CHANNEL
                if (pixelSignMode_ == UNSIGNED_MODE) {
                    rAccumulator += (uint8_t)q.r * (uint8_t)q.a * scaler.r;
                    gAccumulator += (uint8_t)q.g * (uint8_t)q.a * scaler.g;
                    bAccumulator += (uint8_t)q.b * (uint8_t)q.a * scaler.b;
                } else {
                    rAccumulator += (int8_t)q.r * (uint8_t)q.a * scaler.r;
                    gAccumulator += (int8_t)q.g * (uint8_t)q.a * scaler.g;
                    bAccumulator += (int8_t)q.b * (uint8_t)q.a * scaler.b;
                }
        #else
                if (pixelSignMode_ == UNSIGNED_MODE) {
                    rAccumulator += (uint8_t)q.r * scaler.r;
                    gAccumulator += (uint8_t)q.g * scaler.g;
                    bAccumulator += (uint8_t)q.b * scaler.b;
                } else {
                    rAccumulator += (int8_t)q.r * scaler.r;
                    gAccumulator += (int8_t)q.g * scaler.g;
                    bAccumulator += (int8_t)q.b * scaler.b;
                }
        #endif
            }

            // Note: This shift operation will be absolutely necessary when this is implemented
            //       in hardware to avoid the division operation.
        #ifdef USE_ALPHA_CHANNEL
            if (pixelSignMode_ == UNSIGNED_MODE) {
                q.r = CLAMP_UNSIGNED_BYTE(rAccumulator >> (8 + accumulatorShifter_));
                q.g = CLAMP_UNSIGNED_BYTE(gAccumulator >> (8 + accumulatorShifter_));
                q.b = CLAMP_UNSIGNED_BYTE(bAccumulator >> (8 + accumulatorShifter_));
            } else {
                q.r = CLAMP_SIGNED_BYTE(rAccumulator >> (8 + accumulatorShifter_));
                q.g = CLAMP_SIGNED_BYTE(gAccumulator >> (8 + accumulatorShifter_));
                q.b = CLAMP_SIGNED_BYTE(bAccumulator >> (8 + accumulatorShifter_));
            }
        #else
            if (pixelSignMode_ == UNSIGNED_MODE) {
                q.r = CLAMP_UNSIGNED_BYTE(rAccumulator >> accumulatorShifter_);
                q.g = CLAMP_UNSIGNED_BYTE(gAccumulator >> accumulatorShifter_);
                q.b = CLAMP_UNSIGNED_BYTE(bAccumulator >> accumulatorShifter_);
            } else {
                q.r = CLAMP_SIGNED_BYTE(rAccumulator >> accumulatorShifter_);
                q.g = CLAMP_SIGNED_BYTE(gAccumulator >> accumulatorShifter_);
                q.b = CLAMP_SIGNED_BYTE(bAccumulator >> accumulatorShifter_);
            }
        #endif
            q.a = 255;

            if (doCostCalculation) costModel->registerPixelMappingCharge(1);

            frameBuffer_[y * displayWidth_ + x] = q;
        }
    }
}

Pixel* SimpleNddiDisplay::GetFrameBuffer() {
    return GetFrameBuffer(0, 0, displayWidth_, displayHeight_);
}
Pixel* SimpleNddiDisplay::GetFrameBuffer(unsigned int sub_x, unsigned int sub_y, unsigned int sub_w, unsigned int sub_h) {
#ifdef SUPRESS_EXCESS_RENDERING
    if (changed_)
        Render(sub_x, sub_y, sub_w, sub_h);
#endif

    return frameBuffer_;
}
