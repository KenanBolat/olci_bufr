#!/usr/bin/python
# Copyright 2020 EUMETSAT
# Author: Aydin Erturk inspiration from S1 encoder

"""Encode Sentinel 3 OLCI netCDF dataset into BUFR."""

import argparse
import logging
import os
import numpy as np

from netCDF4 import Dataset
from datetime import datetime
from eccodes import (
                     CODES_MISSING_LONG,
                     CODES_MISSING_DOUBLE,
                     codes_bufr_new_from_samples,
                     codes_set_array,
                     codes_set,
                     codes_set_long,
                     codes_set_double_array,
                     codes_write,
                    )

class S3olciBUFR(object):

    """Sentinel 3 OLCI BUFR writer."""

    # Template for output BUFR
    unexpandedDescriptors = [
        1007, 2019, 1096, 25061, 5040, 301011, 301013, 301021, 7025, 5022, 10080, 27080, 
        8003, 8072, 7004, 13093, 8003, 8077, 201131, 202129, 7004, 7004, 202000, 201000, 
        201129, 13095, 201000, 8093, 13095, 8093
        ]

    def __init__(self, infile, outfile):
        """
        Store inputs, extract variables, map them and write to output file.
        """

        self.outfile = outfile

        bufr = self.setup_bufr_metadata()

        logging.info("Processing %s" % infile)
        vals, dims, attrs = self.read_s3olci_netcdf(infile)
        bufr = self.populate_bufr(bufr, vals, dims, attrs)
        self.write_output_bufr(bufr)

    def populate_bufr(self, bufr, vals, dims, attrs):
        """Populate BUFR with values from netCDF."""

        def extract_time(dims, vals):
            """extract time from time stamp"""
            # 

        def extract_metadata(attrs):
            """Interpret attributes for use in BUFR."""
            # Platform
            platform=str(attrs['product_name'][2:5])
            if('A' in platform):
               satelliteID=61
            elif('B' in platform):
               satelliteID=65
            return (platform, satelliteID)

        def set_metadata(bufr, attrs, satelliteID ):
            """Set identifying metadata."""
            # Set metadata
            codes_set(bufr, 'typicalYear', int(attrs['start_time'][2:6]))
            codes_set(bufr, 'typicalMonth', int(attrs['start_time'][7:9]))
            codes_set(bufr, 'typicalDay', int(attrs['start_time'][10:12]))
            codes_set(bufr, 'typicalHour', int(attrs['start_time'][13:15]))
            codes_set(bufr, 'typicalMinute', int(attrs['start_time'][16:18]))
            codes_set(bufr, 'typicalSecond', int(attrs['start_time'][19:21]))
            codes_set(bufr, 'satelliteIdentifier', satelliteID) 
            codes_set(bufr, 'satelliteInstruments', 179)
            codes_set(bufr, 'stationAcquisition', (attrs['institution']))
            #codes_set(bufr, 'softwareVersionNumber', (attrs['source'][11:]))
            codes_set(bufr, 'orbitNumber', int(attrs['absolute_orbit_number']))
            codes_set(bufr, 'year', int(attrs['start_time'][2:6]))
            codes_set(bufr, 'month', int(attrs['start_time'][7:9]))
            codes_set(bufr, 'day', int(attrs['start_time'][10:12]))
            codes_set(bufr, 'hour', int(attrs['start_time'][13:15]))
            codes_set(bufr, 'minute', int(attrs['start_time'][16:18]))
            codes_set(bufr, 'second', int(attrs['start_time'][19:21]))
            return bufr

        def encode_observations(bufr, dims, vals):
            """Encode observations into BUFR."""
            #date = datetime.fromtimestamp(vals['time_stamp'][0]/1000000 + 946681200) # convert milisec to sec / add seconds from year 1900
            for t in xrange(len(dims['rows'])):
                print (t)
                date = datetime.fromtimestamp(vals['time_stamp'][t]/1000000 + 946681200)
                codes_set(bufr, 'year', date.year)
                codes_set(bufr, 'month', date.month)
                codes_set(bufr, 'day', date.day)
                codes_set(bufr, 'hour', date.hour)
                codes_set(bufr, 'minute', date.minute)
                codes_set(bufr, 'second', date.second)
                #codes_set_double_array(bufr, 'longitude(highAccuracy)',vals['longitude'][t])
                #for en in xrange(len(dims['columns'])):
                    #codes_set(bufr, "#%d#latitude(highAccuracy)" %(en+1),vals['latitude'][t][en])
                #codes_set_double_array(bufr, 'solarZenithAngle',vals['SZA'][t])
                #codes_set_double_array(bufr, 'solarAzimuth',vals['SAA'][t])
                #codes_set_double_array(bufr, 'viewingZenithAngle',vals['OZA'][t])
                #codes_set_double_array(bufr, 'viewingAzimuthAngle',vals['OAA'][t])
                codes_set(bufr, 'verticalSignificance(satelliteObservations)',2)
                codes_set(bufr, 'pixel(sType',0)
                codes_set(bufr, 'pressure',CODES_MISSING_DOUBLE)
                codes_set(bufr, 'cloudOpticalThickness',CODES_MISSING_DOUBLE)
                codes_set(bufr, 'verticalSignificance(satelliteObservations)',0)
                #codes_set_array(bufr, '#%d#radiometerSensedSurfaceType',vals['WQFS'][t])
                #codes_set(bufr, '#%d#presure'%(t+1),1)
                #codes_set_double_array(bufr, "#%d#presure"%(p+1),vals['sea_level_pressure'][t])
                #codes_set_double_array(bufr, "#%d#totalColumnWaterVapour"%(t+1),vals['IWV'][t])
                #codes_set(bufr, '#%d#measurementUncertaintySignificance'%(t+1),0)
                #codes_set_double_array(bufr, "#%d#totalColumnWaterVapour"%(t+1),vals['IWV_err'][t])
                #codes_set(bufr, '#%d#measurementUncertaintySignificance'%(t+1),CODES_MISSING_DOUBLE)
            return bufr

        platform, satelliteID = extract_metadata(attrs)

        # Setup structure of BUFR
        codes_set(bufr, 'numberOfSubsets', len(dims['rows']))
        #delayedDescriptorReplication=[len(dims['oswPartitions']),len(dims['oswAngularBinSize'])]
        #delayedDescriptorReplication=len(dims['columns'])
        #codes_set_array(bufr, 'inputDelayedDescriptorReplicationFactor', delayedDescriptorReplication)
        codes_set_array(bufr, 'unexpandedDescriptors', self.unexpandedDescriptors)

        bufr = set_metadata(bufr, attrs, satelliteID)
        bufr = encode_observations(bufr, dims, vals)
        return bufr

    def write_output_bufr(self, bufr):
        """
        Write BUFR constructed in memory to outfile.
        """
        codes_set(bufr, 'pack', True)
        fbufrout = open(self.outfile, 'wb')
        codes_write(bufr, fbufrout)
        logging.info('Created output BUFR file: %s' % self.outfile)
        fbufrout.close()

    def read_s3olci_netcdf(self, s3olci):
        """
        Read input S3 OLCI netCDF file.

        Return dictionaries for variables, dimensions and attributes.
        """
        ds = Dataset(s3olci)
        # read all vaiables, attributes and dimensions in all groups from the nc file and store it in dictionary
        grps = {}
        vals = {}
        dims = {}
        attrs= {}
        for g in ds.groups:
            grps[g] = ds.groups[g]
            for v in ds[g].variables:
                vals[v] = ds[g].variables[v][:]
            for d in ds[g].dimensions:
                dims[d] = ds[g].dimensions[d]
            for a in ds[g].ncattrs():
                attrs[a] = repr(ds[g].getncattr(a))
        # read all attributes from the nc file and store it in dictionary
        #attrs= {}
        #for a in ds.ncattrs():
        #     attrs[a] = repr(ds.getncattr(a))

        return vals, dims, attrs

    def setup_bufr_metadata(self):
        """Create BUFR in memory and insert appropriate metadata."""
        bufr = codes_bufr_new_from_samples('BUFR4')
        if bufr is None:
            logging.info('Not able to codes_bufr_new_from_file')
            raise LookupError("BUFR template not found; check ecCodes paths.")
        codes_set(bufr, 'section2Present', 0)
        codes_set(bufr, 'edition', 4)
        codes_set(bufr, 'masterTableNumber', 0)
        codes_set(bufr, 'bufrHeaderCentre', 254)
        codes_set(bufr, 'bufrHeaderSubCentre', 0)
        codes_set(bufr, 'updateSequenceNumber', 0)
        codes_set(bufr, 'dataCategory', 3)
        codes_set(bufr, 'internationalDataSubCategory', 255)
        codes_set(bufr, 'dataSubCategory', 232)
        codes_set(bufr, 'masterTablesVersionNumber', 34)
        codes_set(bufr, 'localTablesVersionNumber', 0)
        codes_set(bufr, 'observedData', 1)
        codes_set(bufr, 'compressedData', 1)
        return bufr

def set_file_name(inputFile):
    """
    WMO file name convention.
    """
    fileName=os.path.basename(inputFile).split('_',21)
    dateTime=fileName[7].split('T',2)
    satId=fileName[0].upper()
    outFile='W_XX-EUMETSAT-Darmstadt,SURFACE+SATELLITE,'+satId+'+OL+NRT+WV_C_EUMP_'+dateTime[0]+dateTime[1]+'_WV_VV.bin'
    return outFile

def main():
    """
    Instantiate a S3OLCIBUFR object, pass in the netCDF, write it out to file.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="Input netCDF file")
#    parser.add_argument("outfile", help="Output BUFR file")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose mode")
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(format="%(asctime)s: %(message)s", level=logging.INFO)
    outfile = set_file_name(args.infile)
    bufr = S3olciBUFR(args.infile, outfile)

main()