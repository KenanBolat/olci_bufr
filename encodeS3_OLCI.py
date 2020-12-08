#!/usr/bin/python3
# Copyright 2020 EUMETSAT
# Author: Aydin Erturk inspiration from S1 encoder

"""Encode Sentinel 3 OLCI netCDF dataset into BUFR."""

import argparse
import logging
import os
import numpy as np

from netCDF4 import Dataset
from datetime import datetime
# region for debugging
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
# endregion
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

    def __init__(self, infile, outfile, append=False):
        """
        Store inputs, extract variables, map them and write to output file.
        """

        self.outfile = outfile

        bufr = self.setup_bufr_metadata()

        logging.info("Processing %s" % infile)
        vals, dims, attrs = self.read_s3olci_netcdf(infile)
        bufr = self.populate_bufr(bufr, vals, dims, attrs)

    def populate_bufr(self, bufr, vals, dims, attrs):
        """Populate BUFR with values from netCDF."""

        def extract_time(dims, vals):
            """extract time from time stamp"""
            # 

        def extract_metadata(attrs):
            """Interpret attributes for use in BUFR."""
            # Platform
            platform = str(attrs['product_name'][2:5])
            if ('A' in platform):
                satelliteID = 61
            elif ('B' in platform):
                satelliteID = 65
            return (platform, satelliteID)

        def set_metadata(bufr, attrs, satelliteID):
            """Set identifying metadata."""
            # Set metadata
            date_string = attrs['start_time'].replace('\'', '')
            date_value = datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S.%fZ")

            codes_set(bufr, 'typicalYear', date_value.year)
            codes_set(bufr, 'typicalMonth', date_value.month)
            codes_set(bufr, 'typicalDay', date_value.day)
            codes_set(bufr, 'typicalHour', date_value.hour)
            codes_set(bufr, 'typicalMinute', date_value.minute)
            codes_set(bufr, 'typicalSecond', date_value.second)
            codes_set(bufr, 'satelliteIdentifier', satelliteID)
            codes_set(bufr, 'satelliteInstruments', 179)
            codes_set(bufr, 'stationAcquisition', (attrs['institution']))
            # codes_set(bufr, 'softwareVersionNumber', (attrs['source'][11:]))
            codes_set(bufr, 'orbitNumber', int(attrs['absolute_orbit_number']))
            codes_set(bufr, 'year', date_value.year)
            codes_set(bufr, 'month', date_value.month)
            codes_set(bufr, 'day', date_value.day)
            codes_set(bufr, 'hour', date_value.hour)
            codes_set(bufr, 'minute', date_value.minute)
            codes_set(bufr, 'second', date_value.second)
            return bufr

        def encode_observations(bufr, dims, vals):
            """Encode observations into BUFR."""
            SZAintp = np.zeros(vals['longitude'].shape)
            SAAintp = np.zeros(vals['longitude'].shape)
            OZAintp = np.zeros(vals['longitude'].shape)
            OAAintp = np.zeros(vals['longitude'].shape)
            SLPintp = np.zeros(vals['longitude'].shape)
            WQSFintp = np.zeros(vals['longitude'].shape)

            intpFac = ((vals['longitude'].shape[1] - 1) // (vals['SZA'].shape[1] - 1))
            for m in range(vals['SZA'].shape[0] - 14900):
                k = 0
                WQSFintp[m, :] = vals['WQSF'][m, :].astype('int8')
                for i in range(vals['SZA'].shape[1] - 1):
                    for j in range(intpFac):
                        SZAintp[m, k] = vals['SZA'][m, i] + j * ((vals['SZA'][m, i + 1] - vals['SZA'][m, i]) / intpFac)
                        SAAintp[m, k] = vals['SAA'][m, i] + j * ((vals['SAA'][m, i + 1] - vals['SAA'][m, i]) / intpFac)
                        OZAintp[m, k] = vals['OZA'][m, i] + j * ((vals['OZA'][m, i + 1] - vals['OZA'][m, i]) / intpFac)
                        OAAintp[m, k] = vals['OAA'][m, i] + j * ((vals['OAA'][m, i + 1] - vals['OAA'][m, i]) / intpFac)
                        SLPintp[m, k] = vals['sea_level_pressure'][m, i] + j * (
                                (vals['sea_level_pressure'][m, i + 1] - vals['sea_level_pressure'][m, i]) / intpFac)
                        k = k + 1
                        print(m, k, SLPintp[m, k])

            # date = datetime.fromtimestamp(vals['time_stamp'][0]/1000000 + 946681200) # convert milisec to sec / add seconds from year 1900
            scale_factor = 0.299998

            ivw_data = vals['IWV'].filled()
            ivw_data_rectified = np.where((ivw_data < 255) & (ivw_data > 0),
                                          ivw_data * scale_factor,
                                          CODES_MISSING_LONG)

            ivw_err_data = vals['IWV_err'].filled()
            ivw_err_data_rectified = np.where((ivw_err_data < 255) & (ivw_err_data > 0),
                                              ivw_err_data * scale_factor,
                                              CODES_MISSING_LONG)
            SAAintp_rec = np.where((SAAintp < 0 ), SAAintp + 180, SAAintp)
            OAAintp_rec = np.where((OAAintp < 0 ), OAAintp + 180, OAAintp)

            for en, t in enumerate(range(len(dims['rows']))):
                print(t)
                time_stamp_array = datetime.fromtimestamp(vals['time_stamp'][t] / 1000000 + 946681200)
                codes_set_array(bufr, 'longitude(highAccuracy)', vals['longitude'][t, :])
                codes_set_array(bufr, 'latitude(highAccuracy)', vals['latitude'][t, :])
                codes_set_array(bufr, 'solarZenithAngle', SZAintp[t])
                codes_set_array(bufr, 'solarAzimuth', SAAintp_rec[t])
                codes_set_array(bufr, 'viewingZenithAngle', OZAintp[t])
                codes_set_array(bufr, 'viewingAzimuthAngle',OAAintp_rec[t])
                time_stamp_array = datetime.fromtimestamp(vals['time_stamp'][t] / 1000000 + 946681200)
                codes_set(bufr, 'year', time_stamp_array.year)
                codes_set(bufr, 'month', time_stamp_array.month)
                codes_set(bufr, 'day', time_stamp_array.day)
                codes_set(bufr, 'hour', time_stamp_array.hour)
                codes_set(bufr, 'minute', time_stamp_array.minute)
                codes_set(bufr, 'second', time_stamp_array.second)
                codes_set_array(bufr, 'longitude(highAccuracy)', vals['longitude'][t, :])
                codes_set_array(bufr, 'latitude(highAccuracy)', vals['latitude'][t, :])
                codes_set_array(bufr, 'solarZenithAngle', SZAintp[t])
                # codes_set_array(bufr, 'solarAzimuth', SAAintp[t])
                codes_set_array(bufr, 'viewingZenithAngle', OZAintp[t])
                # codes_set_array(bufr, 'viewingAzimuthAngle',OAAintp[t])
                codes_set(bufr, 'verticalSignificance(satelliteObservations)', 2)
                codes_set(bufr, 'pixel(sType', 0)
                codes_set(bufr, 'pressure', CODES_MISSING_DOUBLE)
                codes_set(bufr, 'cloudOpticalThickness', CODES_MISSING_DOUBLE)
                codes_set(bufr, 'verticalSignificance(satelliteObservations)', 0)
                codes_set_array(bufr, 'radiometerSensedSurfaceType', WQSFintp[t, :])
                codes_set(bufr, 'pressure', 1)
                # codes_set(bufr, 'pressure', SLPintp[t,:])
                codes_set_double_array(bufr, "#1#totalColumnWaterVapour", ivw_data_rectified[t, :])
                codes_set_double_array(bufr, "#2#totalColumnWaterVapour", ivw_err_data_rectified[t,: ])
                codes_set(bufr, 'measurementUncertaintySignificance', 0)
                codes_set(bufr, 'measurementUncertaintySignificance', CODES_MISSING_DOUBLE)
                # codes_set(bufr, '#1#pressure', np.ones(vals['longitude'][0,:].shape))
                # codes_set(bufr, '#2#pressure', np.ones(vals['longitude'][0,:].shape)*CODES_MISSING_DOUBLE)
                # codes_set(bufr, '#3#pressure', vals['sea_level_pressure'][t,: ])
                codes_set_array(bufr, 'radiometerSensedSurfaceType', WQSFintp[t, :])

                # Kenan normalde burasi array olacak sunun gibi    codes_set_array(bufr, 'pressure', SLPintp[t,:])

                codes_set(bufr, 'pressure', SLPintp[t, 0] * 100)

                # #codes_set_double_array(bufr, "#%d#pressure"%(p+1),vals['sea_level_pressure'][t])

                # for m in range(len(dims['columns'])):
                #     date = datetime.fromtimestamp(vals['time_stamp'][t]/1000000 + 946681200)
                #     codes_set(bufr, 'year', date.year)
                #     codes_set(bufr, 'month', date.month)
                #     codes_set(bufr, 'day', date.day)
                #     codes_set(bufr, 'hour', date.hour)
                #     codes_set(bufr, 'minute', date.minute)
                #     codes_set(bufr, 'second', date.second)
                #     # codes_set(bufr, 'longitude(highAccuracy)',vals['longitude'][t][m])
                #     # codes_set(bufr, 'latitude(highAccuracy)',vals['latitude'][t][m])
                #     # print (vals['longitude'][t][m])
                # #for en in xrange(len(dims['columns'])):
                #     #codes_set(bufr, "#%d#latitude(highAccuracy)" %(en+1),vals['latitude'][t][en])
                # #codes_set_double_array(bufr, 'solarZenithAngle',vals['SZA'][t])
                # #codes_set_double_array(bufr, 'solarAzimuth',vals['SAA'][t])
                # #codes_set_double_array(bufr, 'viewingZenithAngle',vals['OZA'][t])
                # #codes_set_double_array(bufr, 'viewingAzimuthAngle',vals['OAA'][t])
                #     codes_set(bufr, 'verticalSignificance(satelliteObservations)',2)
                #     codes_set(bufr, 'pixel(sType',0)
                #     codes_set(bufr, 'pressure',CODES_MISSING_DOUBLE)
                #     codes_set(bufr, 'cloudOpticalThickness',CODES_MISSING_DOUBLE)
                #     codes_set(bufr, 'verticalSignificance(satelliteObservations)',0)
                # #codes_set_array(bufr, '#%d#radiometerSensedSurfaceType',vals['WQFS'][t])
                # #codes_set(bufr, '#%d#presure'%(t+1),1)
                # #codes_set_double_array(bufr, "#%d#presure"%(p+1),vals['sea_level_pressure'][t])
                # #codes_set_double_array(bufr, "#%d#totalColumnWaterVapour"%(t+1),vals['IWV'][t])
                # #codes_set(bufr, '#%d#measurementUncertaintySignificance'%(t+1),0)
                # #codes_set_double_array(bufr, "#%d#totalColumnWaterVapour"%(t+1),vals['IWV_err'][t])
                # #codes_set(bufr, '#%d#measurementUncertaintySignificance'%(t+1),CODES_MISSING_DOUBLE)
                if t == 0:
                    append = False
                else:
                    append = True
                self.write_output_bufr(bufr, append)

        platform, satelliteID = extract_metadata(attrs)

        # Setup structure of BUFR
        codes_set(bufr, 'numberOfSubsets', len(dims['columns']))
        # delayedDescriptorReplication=[len(dims['oswPartitions']),len(dims['oswAngularBinSize'])]
        # delayedDescriptorReplication=len(dims['columns'])
        # codes_set_array(bufr, 'inputDelayedDescriptorReplicationFactor', delayedDescriptorReplication)
        codes_set_array(bufr, 'unexpandedDescriptors', self.unexpandedDescriptors)

        bufr = set_metadata(bufr, attrs, satelliteID)
        bufr = encode_observations(bufr, dims, vals)
        return bufr

    def write_output_bufr(self, bufr, append=False):
        """
        Write BUFR constructed in memory to outfile.
        """
        codes_set(bufr, 'pack', True)
        if not append:
            fbufrout = open(self.outfile, 'wb')
        else:
            fbufrout = open(self.outfile, 'ab')
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
        attrs = {}
        for g in ds.groups:
            grps[g] = ds.groups[g]
            for v in ds[g].variables:
                vals[v] = ds[g].variables[v][:]
            for d in ds[g].dimensions:
                dims[d] = ds[g].dimensions[d]
            for a in ds[g].ncattrs():
                attrs[a] = repr(ds[g].getncattr(a))
        # read all attributes from the nc file and store it in dictionary
        # attrs= {}
        # for a in ds.ncattrs():
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
    fileName = os.path.basename(inputFile).split('_', 21)
    dateTime = fileName[7].split('T', 2)
    satId = fileName[0].upper()
    outFile = 'W_XX-EUMETSAT-Darmstadt,SURFACE+SATELLITE,' + satId + '+OL+NRT+WV_C_EUMP_' + dateTime[0] + dateTime[
        1] + '_WV_VV.bin'
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
