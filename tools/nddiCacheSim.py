#!/usr/bin/env python

import sys, argparse
import urllib
from urllib import urlopen
import json
import ijson.backends.yajl2_cffi as ijson
from cachesim import CacheSimulator, Cache, MainMemory

parser = argparse.ArgumentParser(description='Process json from the NDDI cost model and simulate the cache accesses.')
parser.add_argument('file', help='path to the json log file from the NDDI cost model')
parser.add_argument('--strides', help='provide a stride order other than the default 0,1,... (e.g. --strides 2,0,1)')
parser.add_argument('--limit', help='set the limit to the number of charges parsed (useful for debugging)', type=int)
parser.add_argument('--verbose', action='store_const', const=True, default=False,
                    help='prints additional information while processing the json')
parser.add_argument('--ijson', action='store_const', const=True, default=False,
                    help='Uses ijson instead of the json package')

args = parser.parse_args()
if args.verbose:
    print "Arguments: ", args

print "Parsing json from file: ", args.file
fp = urlopen(args.file)
if args.ijson:
    objects = ijson.items(fp, 'config')
    config = objects.next()
else:
    data = json.load(fp)
    config = data['config']
fp.close()

ivSize = config['inputVectorSize']
fvDimensions = config['fvDimensions']
fvStrideOrder = []
if (args.strides == None):
    for idx, val in enumerate(fvDimensions):
        fvStrideOrder.append(idx)
else:
    strides = args.strides.split(',')
    for stride in strides:
        fvStrideOrder.append(int(stride))
bpp = config['bytePerPixel']
bpiv = config['bytePerIvValue']
bpc = config['bytePerCoefficient']
bps = config['bytePerScaler']

print "Input Vector Size: ", ivSize
print "Frame Volume Dimensions: ", fvDimensions
print "Coefficient Matrix Size: ", ivSize, "x", len(fvDimensions)
print "Frame Volume Stride Order: ", fvStrideOrder
print "Bytes per Pixel: ", bpp
print "Bytes per Input Vector value: ", bpiv
print "Bytes per Coefficient: ", bpc
print "Bytes per Scaler: ", bps


mem = MainMemory()
l3 = Cache("L3", 20480, 16, 64, "LRU")  # 20MB: 20480 sets, 16-ways with cacheline size of 64 bytes
mem.load_to(l3)
mem.store_from(l3)
l2 = Cache("L2", 512, 8, 64, "LRU", store_to=l3, load_from=l3)  # 256KB
l1 = Cache("L1", 64, 8, 64, "LRU", store_to=l2, load_from=l2)  # 32KB
cs = CacheSimulator(l1, mem)

def fvAccessRow(tuple, length, access):
    strideMultiplier = bpp
    mem = 0
    bytes = length * bpp
    for strideOrder in fvStrideOrder:
        mem += tuple[strideOrder] * strideMultiplier
        strideMultiplier *= fvDimensions[strideOrder]
    if args.verbose:
        print "Addr: ", mem, " Byte Count: ", bytes
    if access == "READ_ACCESS":
        cs.load(mem, bytes)
    elif access == "WRITE_ACCESS":
        cs.store(mem, bytes)
    else:
        sys.exit('ERROR Unknown access type. Should be READ_ACCESS or WRITE_ACCESS only.')

charges = []
if args.ijson:
    fp = urlopen(args.file)
    charges = ijson.items(fp, 'charges.item')
else:
    charges = data['charges']

count = 0
for charge in charges:
    count += 1
    if args.limit and count == args.limit:
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
                    fvAccessRow(current, end[strideOrder] - start[strideOrder] + 1, access)
                    carry = 1
                else:
                    current[strideOrder] += carry
                    if current[strideOrder] > end[strideOrder]:
                        current[strideOrder] = start[strideOrder]
                        carry = 1
                    else:
                        carry = 0
    else:
        sys.exit('ERROR Can only handle Frame Volume Charges at this time.')

print count
if args.ijson:
    fp.close()

cs.force_write_back()
cs.print_stats()
