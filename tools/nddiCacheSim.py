#!/usr/bin/env python

import sys, json, argparse
from cachesim import CacheSimulator, Cache, MainMemory

parser = argparse.ArgumentParser(description='Process json from the NDDI cost model and simulate the cache accesses.')
parser.add_argument('file', help='path to the json log file from the NDDI cost model')
parser.add_argument('--strides', help='provide a stride order other than the default 0,1,... (e.g. --strides 2,0,1)')
parser.add_argument('--limit', help='set the limit to the number of charges parsed (useful for debugging)', type=int)
parser.add_argument('--verbose', action='store_const', const=True, default=False,
                    help='prints additional information while processing the json')

args = parser.parse_args()
if args.verbose:
    print "Arguments: ", args

print "Parsing json from file: ", args.file
fp = open(args.file)
data = json.load(fp)
fp.close()

ivSize = data["inputVectorSize"]
fvDimensions = data["fvDimensions"]
fvStrideOrder = []
if (args.strides == None):
    for idx, val in enumerate(fvDimensions):
        fvStrideOrder.append(idx)
else:
    strides = args.strides.split(',')
    for stride in strides:
        fvStrideOrder.append(int(stride))
print "Input Vector Size: ", ivSize
print "Frame Volume Dimensions: ", fvDimensions
print "Coefficient Matrix Size: ", ivSize, "x", len(fvDimensions)
print "Frame Volume Stride Order: ", fvStrideOrder
print "Bytes per Pixel: ", data["bytePerPixel"]
print "Bytes per Input Vector value: ", data["bytePerIvValue"]
print "Bytes per Coefficient: ", data["bytePerCoefficient"]
print "Bytes per Scaler: ", data["bytePerScaler"]


mem = MainMemory()
l3 = Cache("L3", 20480, 16, 64, "LRU")  # 20MB: 20480 sets, 16-ways with cacheline size of 64 bytes
mem.load_to(l3)
mem.store_from(l3)
l2 = Cache("L2", 512, 8, 64, "LRU", store_to=l3, load_from=l3)  # 256KB
l1 = Cache("L1", 64, 8, 64, "LRU", store_to=l2, load_from=l2)  # 32KB
cs = CacheSimulator(l1, mem)

def fvToMem(a):
    strideMultiplier = data["bytePerPixel"]
    m = 0
    for strideOrder in fvStrideOrder:
        m += a[strideOrder] * strideMultiplier
        strideMultiplier *= fvDimensions[strideOrder]
    return m

count = args.limit
for charge in data["charges"]:
    count = count - 1
    if count < 0:
        break
    if "frameVolumeCharge" in charge:
        if args.verbose:
            print charge["frameVolumeCharge"]
        
        access = charge["frameVolumeCharge"]["access"]
        start = charge["frameVolumeCharge"]["start"]
        end = charge["frameVolumeCharge"]["end"]
        
        carry = 0
        current = list(start)
        while carry == 0:
            for idx, strideOrder in enumerate(fvStrideOrder):
                if idx == 0:
                    mem = fvToMem(current)
                    bytes = (end[strideOrder] - start[strideOrder] + 1) * data["bytePerPixel"]
                    carry = 1
                    if args.verbose:
                        print "Addr: ", mem, " Byte Count: ", bytes
                    if access == "READ_ACCESS":
                        cs.load(mem, length = bytes)
                    elif access == "WRITE_ACCESS":
                        cs.store(mem, length = bytes)
                    else:
                        sys.exit('ERROR Unknown access type. Should be READ_ACCESS or WRITE_ACCESS only.')
                else:
                    current[strideOrder] += carry
                    if current[strideOrder] > end[strideOrder]:
                        current[strideOrder] = start[strideOrder]
                        carry = 1
                    else:
                        carry = 0
    else:
        sys.exit('ERROR Can only handle Frame Volume Charges at this time.')

cs.force_write_back()
cs.print_stats()
