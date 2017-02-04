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
        : size_(size), values_(NULL), costModel_(costModel) {

            if (!costModel_->isHeadless()) {
                values_ = (int *)malloc(sizeof(int) * size);
                memset(values_, 0x00, sizeof(int) * size);
            }
        }

        ~InputVector() {

            if (values_) {
                free(values_);
            }
        }

        unsigned int getSize() {

            return size_;
        }

        void UpdateInputVector(int* input) {

            if (!costModel_->isHeadless()) {
                for (int i = 0; (i + 2) < size_; i++) {
                    setValue(i+2, input[i]);
                }
            } else {
                costModel_->registerBulkMemoryCharge(INPUT_VECTOR_COMPONENT,
                                                     size_,
                                                     WRITE_ACCESS,
                                                     NULL,
                                                     size_ * BYTES_PER_IV_VALUE,
                                                     0);
            }
        }

        void setValue(unsigned int location, int value) {

            assert(location < size_);
            assert(!costModel_->isHeadless());

            costModel_->registerMemoryCharge(INPUT_VECTOR_COMPONENT, WRITE_ACCESS, values_ + location, BYTES_PER_IV_VALUE, 0);

            values_[location] = value;
        }

        int getValue(unsigned int location) {

            assert(location > 1);
            assert(location < size_);
            assert(!costModel_->isHeadless());

            costModel_->registerMemoryCharge(INPUT_VECTOR_COMPONENT, READ_ACCESS, values_ + location, BYTES_PER_IV_VALUE, 0);

            return values_[location];
        }

        int * data() {

            return values_;
        }

    };
}

#endif
