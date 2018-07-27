//
//  InputVector.h
//  pixelbridge
//
//  Created by Dave Estes on 2/29/12.
//  Copyright 2012 Dave Estes. All rights reserved.
//

#ifndef pixelbridge_InputVector_h
#define pixelbridge_InputVector_h

#include <string.h>
#include <cstdlib>
#include <cassert>

#include "NDimensionalDisplayInterface.h"

using namespace std;

namespace nddi {

    class InputVector {

    protected:
        CostModel     * costModel_;
        size_t          size_;
        int           * values_;

    public:

        InputVector(CostModel* costModel,
                    unsigned int size)
        : size_(0), values_(NULL), costModel_(costModel) {

            if (!costModel_->isHeadless()) {
                values_ = (int *)malloc(sizeof(int) * size);
                memset(values_, 0x00, sizeof(int) * size);
            }
            size_ = size;
        }

        ~InputVector() {

            if (values_) {
                free(values_);
            }
        }

        unsigned int getSize() {

            return size_;
        }

        void UpdateInputVector(vector<int> &input) {

            assert(input.size() + 2 == size_);

            if (!costModel_->isHeadless()) {
                for (int i = 0; (i < input.size()) && ((i + 2) < size_); i++) {
                    setValue(i+2, input[i]);
                }
            }
            costModel_->registerInputVectorMemoryCharge(WRITE_ACCESS, 2, input.size() - 1);
        }

        void setValue(unsigned int location, int value) {

            assert(location < size_);
            assert(!costModel_->isHeadless());

            values_[location] = value;
        }

        int getValue(unsigned int location) {

            assert(location > 1);
            assert(location < size_);
            assert(!costModel_->isHeadless());

            return values_[location];
        }

        int * data() {

            return values_;
        }

    };
}

#endif
