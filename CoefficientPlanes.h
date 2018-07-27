//
//  CoefficientPlanes.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_CoefficientPlanes_h
#define pixelbridge_CoefficientPlanes_h

#include <algorithm>
#include <cstdlib>
#include <cassert>

#include "NDimensionalDisplayInterface.h"
#include "CoefficientMatrix.h"

/**
 * These helper macros will calculate the offset into scaler memory or coefficient memory. This is helpful
 * for tweaking the layout when doing memory profiling experiments.
 */
#define MB_FACTOR         (fixed8x8Macroblocks_ ? 8 : 1)
#define MB_SCALE(v)       (fixed8x8Macroblocks_ ? v >> 3 : v)
#define MB_SCALE_UP(v)    (fixed8x8Macroblocks_ ? (v + 7) >> 3 : v)
#define NUM_COEFF_PLANES  (useSingleCoefficientPlane_ ? 1 : numPlanes_)

// R x C x P
#define SC_OFF(x, y, p, c)   ((MB_SCALE(y) * MB_SCALE_UP(width_) * numPlanes_ + MB_SCALE(x) * numPlanes_ + p) * 3 + c)
#define CP_OFF(x, y, p)   ((MB_SCALE(y) * MB_SCALE_UP(width_) * NUM_COEFF_PLANES + MB_SCALE(x) * NUM_COEFF_PLANES + p) * matrixWidth_ * matrixHeight_)

// P x R x C
//#define SC_OFF(x, y, p, c)   ((p * MB_SCALE_UP(height_) * MB_SCALE_UP(width_) + MB_SCALE(y) * MB_SCALE_UP(width_) + MB_SCALE(x)) * 3 + c)
//#define CP_OFF(x, y, p)   ((p * MB_SCALE_UP(height_) * MB_SCALE_UP(width_) + MB_SCALE(y) * MB_SCALE_UP(width_) + MB_SCALE(x)) * matrixWidth_ * matrixHeight_)

using namespace std;

namespace nddi {

    class CoefficientPlanes {

    protected:
        CostModel           * costModel_;
        size_t                width_, height_, numPlanes_, matrixWidth_, matrixHeight_;
        bool                  fixed8x8Macroblocks_;
        bool                  useSingleCoefficientPlane_;
        CoefficientMatrix   * coefficientMatrix_;
        Coeff               * coefficients_;
        int16_t             * scalers_;

    public:

        CoefficientPlanes() {
        }

        CoefficientPlanes(CostModel* costModel,
                         unsigned int displayWidth, unsigned int displayHeight,
                         unsigned int numPlanes,
                         unsigned int matrixWidth, unsigned int matrixHeight,
                         bool fixed8x8Macroblocks = false, bool useSingleCoefficientPlane = false)
        : costModel_(costModel),
          width_(displayWidth), height_(displayHeight),
          numPlanes_(numPlanes),
          matrixWidth_(matrixWidth), matrixHeight_(matrixHeight),
          fixed8x8Macroblocks_(fixed8x8Macroblocks),
          useSingleCoefficientPlane_(useSingleCoefficientPlane) {

            // Create the common CoefficientMatrix
            coefficientMatrix_ =  new CoefficientMatrix(costModel_, matrixWidth, matrixHeight);

            // Alloc the actual coefficients and scalers
            if (!costModel_->isHeadless()) {
#ifdef DEBUG
                coefficients_ = (Coeff *)calloc(1, CoefficientMatrix::memoryRequired(matrixWidth, matrixHeight) * MB_SCALE_UP(displayWidth) * MB_SCALE_UP(displayHeight) * (useSingleCoefficientPlane_ ? 1 : numPlanes_));
                scalers_ = (int16_t *)calloc(1, sizeof(int16_t) * 3 * MB_SCALE_UP(displayWidth) * MB_SCALE_UP(displayHeight) * numPlanes_);
#else
                coefficients_ = (Coeff *)malloc(CoefficientMatrix::memoryRequired(matrixWidth, matrixHeight) * MB_SCALE_UP(displayWidth) * MB_SCALE_UP(displayHeight) * (useSingleCoefficientPlane_ ? 1 : numPlanes_));
                scalers_ = (int16_t *)malloc(sizeof(int16_t) * 3 * MB_SCALE_UP(displayWidth) * MB_SCALE_UP(displayHeight) * numPlanes_);
#endif
            }
        }

        ~CoefficientPlanes() {

            if (coefficientMatrix_) delete(coefficientMatrix_);
            if (coefficients_) free(coefficients_);
            if (scalers_) free(scalers_);
        }

        unsigned int getWidth() {

            return width_;
        }

        unsigned int getHeight() {

            return height_;
        }

        Coeff CheckSpecialCoefficient(Coeff c, unsigned int p) {
            Coeff ret = c;
            switch (c) {
            case COEFFICIENT_MATRIX_X:
                assert(false);
                break;
            case COEFFICIENT_MATRIX_Y:
                assert(false);
                break;
            case COEFFICIENT_MATRIX_P:
                ret = p;
                break;
            default:
                break;
            }

            return ret;
        }

        void PutCoefficientMatrix(vector< vector<int> > &coefficientMatrix, vector<unsigned int> &location) {

            assert(location.size() == 3);
            assert(location[0] < width_);
            assert(location[1] < height_);
            assert(location[2] < numPlanes_);
            assert(coefficientMatrix.size() == matrixWidth_);
            assert(coefficientMatrix[0].size() == matrixHeight_);
            assert(!fixed8x8Macroblocks_ || (!(location[0] % 8) && !(location[1] % 8)));
            assert(!useSingleCoefficientPlane_ || location[2] == 0);

            if (!costModel_->isHeadless()) {
                // Examine each coefficient in the coefficient matrix vector and use it unless it's a COEFFICIENT_UNCHANGED
                for (int col = 0; col < matrixHeight_; col++) {
                    for (int row = 0; row < matrixWidth_; row++) {
                        if (coefficientMatrix[row][col] != COEFFICIENT_UNCHANGED) {
    #ifdef NARROW_DATA_STORES
                            assert(coefficientMatrix[row][col] >= SHRT_MIN && coefficientMatrix[row][col] <= SHRT_MAX);
    #endif
                            Coeff * cm = dataCoefficient(location[0], location[1], location[2]);
                            cm[col * matrixWidth_ + row] = coefficientMatrix[row][col];
                        }
                    }
                }
            }
            costModel_->registerCoefficientMatrixMemoryCharge(WRITE_ACCESS, location, location, coefficientMatrix);
        }

        Coeff getCoefficient(vector<unsigned int> &location, int row, int col) {
            assert(location.size() == 3);
            assert(location[0] < width_);
            assert(location[1] < height_);
            assert(location[2] < numPlanes_);
            assert(row < matrixWidth_);
            assert(col < matrixHeight_);

            Coeff *cm = dataCoefficient(location[0], location[1], useSingleCoefficientPlane_ ? 0 : location[2]);

            return CheckSpecialCoefficient(cm[col * matrixWidth_ + row], location[2]);
        }

        void FillCoefficientMatrix(vector< vector<int> > &coefficientMatrix,
                                   vector<unsigned int> &start,
                                   vector<unsigned int> &end) {

            assert(start.size() == 3);
            assert(start[0] < width_);
            assert(start[1] < height_);
            assert(start[2] < numPlanes_);
            assert(end.size() == 3);
            assert(end[0] < width_);
            assert(end[1] < height_);
            assert(end[2] < numPlanes_);
            assert(coefficientMatrix.size() == matrixWidth_);
            assert(coefficientMatrix[0].size() == matrixHeight_);
            // Note: Unlike the other routines, this one doesn't assert for planes >0 when using useSingleCoefficientPlane_.
            //       Instead it handles those planes for purposes of cost model calculations without actually updating memory.
            // assert(!useSingleCoefficientPlane_ || (start[2] == 0 && end[2] == 0));

            vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Update coefficient matrix in coefficient plane at position
                if (!costModel_->isHeadless()) {
                    // Examine each coefficient in the coefficient matrix vector and use it unless it's a COEFFICIENT_UNCHANGED
                    for (int col = 0; col < matrixHeight_; col++) {
                        for (int row = 0; row < matrixWidth_; row++) {
                            if (coefficientMatrix[row][col] != COEFFICIENT_UNCHANGED) {
        #ifdef NARROW_DATA_STORES
                                assert(coefficientMatrix[row][col] >= SHRT_MIN && coefficientMatrix[row][col] <= SHRT_MAX);
        #endif
                                Coeff * cm = !useSingleCoefficientPlane_ ?
                                        dataCoefficient(position[0], position[1], position[2]) :
                                        dataCoefficient(position[0], position[1], 0);
                                // Only physically update the coefficients on plane zero when using useSingleCoefficientPlane_
                                if (!useSingleCoefficientPlane_ || (position[2] == 0)) {
                                    cm[col * matrixWidth_ + row] = coefficientMatrix[row][col];
                                }
                            }
                        }
                    }
                }

                // Move to the next position
                position[0] += MB_FACTOR;
                if (position[0] > end[0]) {
                    position[0] = start[0];
                    position[1] += MB_FACTOR;
                    if (position[1] > end[1]) {
                        position[1] = start[1];
                        position[2]++;
                        if (position[2] > end[2]) {
                            fillFinished = true;
                        }
                    }
                }
            } while (!fillFinished);

            costModel_->registerCoefficientMatrixMemoryCharge(WRITE_ACCESS, start, end, coefficientMatrix);
        }

        void FillCoefficient(int coefficient,
                             int row, int col,
                             vector<unsigned int> &start,
                             vector<unsigned int> &end) {

            assert(start.size() == 3);
            assert(start[0] < width_);
            assert(start[1] < height_);
            assert(start[2] < numPlanes_);
            assert(end.size() == 3);
            assert(end[0] < width_);
            assert(end[1] < height_);
            assert(end[2] < numPlanes_);
            assert(col < matrixWidth_);
            assert(row < matrixHeight_);
#ifdef NARROW_DATA_STORES
            assert(coefficient >= SHRT_MIN && coefficient <= SHRT_MAX);
#endif
            // Note: Unlike the other routines, this one doesn't assert for planes >0 when using useSingleCoefficientPlane_.
            //       Instead it handles those planes for purposes of cost model calculations without actually updating memory.
            // assert(!useSingleCoefficientPlane_ || (start[2] == 0 && end[2] == 0));

            vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set coefficient in the coefficient matrix at this position in the coefficient plane
                if (!costModel_->isHeadless()) {
                    Coeff * cm = !useSingleCoefficientPlane_ ?
                            dataCoefficient(position[0], position[1], position[2]) :
                            dataCoefficient(position[0], position[1], 0);
                    // Only physically update the coefficients on plane zero when using useSingleCoefficientPlane_
                    if (!useSingleCoefficientPlane_ || (position[2] == 0)) {
                        cm[row * matrixWidth_ + col] = coefficient;
                    }
                }

                // Move to the next position
                position[0] += MB_FACTOR;
                if (position[0] > end[0]) {
                    position[0] = start[0];
                    position[1] += MB_FACTOR;
                    if (position[1] > end[1]) {
                        position[1] = start[1];
                        position[2]++;
                        if (position[2] > end[2]) {
                            fillFinished = true;
                        }
                    }
                }
            } while (!fillFinished);

            costModel_->registerCoefficientMemoryCharge(WRITE_ACCESS, start, end, row, col);
        }

        void FillScaler(Scaler scaler,
                        vector<unsigned int> &start,
                        vector<unsigned int> &end) {

            assert(start.size() == 3);
            assert(start[0] < width_);
            assert(start[1] < height_);
            assert(start[2] < numPlanes_);
            assert(end.size() == 3);
            assert(end[0] < width_);
            assert(end[1] < height_);
            assert(end[2] < numPlanes_);

            vector<unsigned int> position = start;
            bool fillFinished = false;

            // Move from start to end, filling in each location with the provided pixel
            do {
                // Set scaler at this position in the coefficient plane
                if (!costModel_->isHeadless()) {
                    putScaler(position[0], position[1], position[2], scaler);
                }

                // Move to the next position
                position[0] += MB_FACTOR;
                if (position[0] > end[0]) {
                    position[0] = start[0];
                    position[1] += MB_FACTOR;
                    if (position[1] > end[1]) {
                        position[1] = start[1];
                        position[2]++;
                        if (position[2] > end[2]) {
                            fillFinished = true;
                        }
                    }
                }
            } while (!fillFinished);

            costModel_->registerScalerMemoryCharge(WRITE_ACCESS, start, end);
        }

        void FillScalerStack(vector<uint64_t> &scalers,
                             vector<unsigned int> &start,
                             vector<unsigned int> &size) {

        }

        CoefficientMatrix* getCoefficientMatrix() {

            return coefficientMatrix_;
        }

        void putScaler(unsigned int x, unsigned int y, unsigned int p, Scaler scaler) {

            assert(x < width_);
            assert(y < height_);
            assert(!costModel_->isHeadless());
            assert(!fixed8x8Macroblocks_ || (!(x % 8) && !(y % 8)));

            scalers_[SC_OFF(x, y, p, 0)] = scaler.r;
            scalers_[SC_OFF(x, y, p, 1)] = scaler.g;
            scalers_[SC_OFF(x, y, p, 2)] = scaler.b;
        }

        Scaler getScaler(unsigned int x, unsigned int y, unsigned int p) {

            assert(x < width_);
            assert(y < height_);
            assert(p < numPlanes_);
            assert(!costModel_->isHeadless());

            Scaler s;
            s.packed = 0;
            s.r = scalers_[SC_OFF(x, y, p, 0)];
            s.g = scalers_[SC_OFF(x, y, p, 1)];
            s.b = scalers_[SC_OFF(x, y, p, 2)];

            return s;
        }

        int16_t * dataScaler(size_t x, size_t y, size_t p) {
            assert(!costModel_->isHeadless());
            return &scalers_[SC_OFF(x, y, p, 0)];
        }

        Coeff * dataCoefficient(size_t x, size_t y, size_t p) {
            assert(!costModel_->isHeadless());
            assert(!useSingleCoefficientPlane_ || p == 0);
            return &coefficients_[CP_OFF(x, y, p)];
        }

    };
}

#endif
