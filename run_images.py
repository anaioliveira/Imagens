##################################################################
#
#     Developed by: Ana Isabel Oliveira
#     Project: Water4Ever
#     Date: MARETEC IST, 10/05/2019
#     Description: Subtract the results of MOHID to satellite images
#     Outputs: Raster files with the results of the difference
#
##################################################################


#!/usr/bin/python
# -*- coding: utf-8 -*-

# Imports
import os
import sys
import h5py
import numpy
import subprocess
import gdal
from scipy.interpolate import RegularGridInterpolator
import math

def construct_ascii_grid_file(lat, lon, values):

    origin_x = lon.min()
    origin_y = lat.min()
    
    cell_size_lat = round(lat[1]-lat[0], 6)
    cell_size_lon = round(lon[1]-lon[0], 6)
    
    if cell_size_lat != cell_size_lon:
        print ("Cell size in xx axis is different of cell size in yy axis!")
        sys.exit()
    
    ncols = len(lon) - 1
    nrows = len(lat) - 1
    
    values_mohid = numpy.zeros((ncols, nrows))

    file_in = "griddata.asc"
    
    with open("griddata.asc","w") as fin_asc:
        fin_asc.write('ncols ' + str(ncols) + '\n')
        fin_asc.write('nrows ' + str(nrows) + '\n')
        fin_asc.write('xllcorner ' + str(origin_x) + '\n')
        fin_asc.write('yllcorner ' + str(origin_y) + '\n')
        fin_asc.write('cellsize ' + str(cell_size_lat) + '\n')
        fin_asc.write('nodata_value ' + str(values[0][0]) + '\n')
        for r in range(nrows):
            line_to_write = ''
            for c in range(ncols):
                line_to_write = line_to_write + ' ' + str(values[c][nrows-r-1])
            fin_asc.write(line_to_write + '\n')

    fin_asc.close()

    return file_in

if __name__ == '__main__':

    satellite_images = ['A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2015-06-10_end_2015-06-10/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332015161_20150610_111345_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2015-06-26_end_2015-06-26/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332015177_20150626_111350_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2015-07-12_end_2015-07-12/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332015193_20150712_111402_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2015-07-28_end_2015-07-28/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332015209_20150728_111407_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2017-07-01_end_2017-07-01/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332017182_20170701_111417_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2017-07-17_end_2017-07-17/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332017198_20170717_111421_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2017-08-02_end_2017-08-02/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332017214_20170802_111429_LAI.tif',\
                        'A:/Ana/Aplica/Scripts/Imagens/SatelliteImages/mosaic_start_2017-08-18_end_2017-08-18/0476f029-c6b3-40ed-8752-98234deeb2da__landsat8_2040332017230_20170818_111434_LAI.tif'\
                        ]
                        
    hdf_files = ['A:/Ana/Aplica/Scripts/Imagens/HDF/20150610/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20150626/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20150712/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20150728/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20170701/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20170717/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20170802/Vegetation_4.hdf5',\
                 'A:/Ana/Aplica/Scripts/Imagens/HDF/20170818/Vegetation_4.hdf5'\
                 ]
    
    property = 'leaf area index'
    
    instants = ['2015-06-10',\
                '2015-06-26',\
                '2015-07-12',\
                '2015-07-28',\
                '2017-07-01',\
                '2017-07-17',\
                '2017-08-02',\
                '2017-08-18'
                ]
    
    if len(satellite_images) != len(hdf_files):
        print ('Please check if you have the same number of images and HDF5 files to compare.')
        
    if len(hdf_files) != len(instants):
        print ('Please check if you have the same number of instants and HDF5 files.')
    
    #Create the output file
    fi= open("Statistics.dat","w+")
    fi.write("Date R2 RMSE PBias\n")                       
    f = 0
    for instant in instants:
        print ("Working on instant: " + instant)
        hdf5_in_h5py = h5py.File(hdf_files[f], 'r')
        time_keys = list(hdf5_in_h5py['Time'].keys())
        subdatasets = list(hdf5_in_h5py['Results'][property].keys())
        
        instant_aux = instant.split('-')
        year = float(instant_aux[0])
        month = float(instant_aux[1])
        day = float(instant_aux[2])
        
        for t in time_keys:
            t_ = hdf5_in_h5py['Time'][t].value
            if year in t_ and month in t_ and day in t_:
                time_instant = t.split('_')[1]
                break
            else:
                pass
        
        # Get the right dataset
        for dataset in subdatasets:
            if property in dataset and time_instant in dataset:
                values = numpy.array(hdf5_in_h5py['Results'][property][dataset].value)
        
        lat_mohid = hdf5_in_h5py['Grid']['Latitude'].value[0,:]
        lon_mohid = hdf5_in_h5py['Grid']['Longitude'].value[:,0]
        
        ascii_from_hdf= construct_ascii_grid_file(lat_mohid, lon_mohid, values)
        
        tif_from_hdf = ascii_from_hdf.split('.')[0]+'.tif'
        
        os.system('gdal_translate -of "GTiff" -a_srs EPSG:4326 ' + ascii_from_hdf + " " + tif_from_hdf)
        
        #--------------------Satellite images---------------------
        satellite_image = satellite_images[f]
        file_out = "difference_" + instant + ".tif"
        
        satellite_image_wgs = satellite_image.split('.')[0] + '_wgs84.tif'
        
        gdal.Warp(satellite_image_wgs,satellite_image,dstSRS='EPSG:4326')
        
        # Get satellite image data
        ds1_src = gdal.Open(satellite_image_wgs)
        colum_numbers,row_numbers=ds1_src.RasterXSize,ds1_src.RasterYSize
        
        # GDAL affine transform parameters, According to gdal documentation xoff/yoff are image left corner, a/e are pixel wight/height and b/d is rotation and is zero if image is north up.
        xoff, a, b, yoff, d, e = ds1_src.GetGeoTransform()
        
        #Compute latitude and longitude of corners
        lon_ds1_ = []
        lat_ds1_ = []
        for c in range(row_numbers+1):
            lat_ds1_.append(yoff + e * c)
        lat_ds1_.sort()

        for r in range(colum_numbers+1):
            lon_ds1_.append(xoff + a * r)
        lon_ds1_.sort()
            
        #Copute latitude and logitude of cell centers
        lon_ds1 = []
        lat_ds1 = []
        for k in range(1,len(lat_ds1_)):
            lat_ds1.append((lat_ds1_[k]+lat_ds1_[k-1])/2)

        for m in range(1,len(lon_ds1_)):
            lon_ds1.append((lon_ds1_[m]+lon_ds1_[m-1])/2)
            
        # Dealing with the band values
        ds1_bnd_ = ds1_src.GetRasterBand(1).ReadAsArray()/1000
        ds1_bnd_ [ds1_bnd_ > 50] = -99
        ds1_bnd = numpy.flip(ds1_bnd_, axis=0)

        interpolated_satellite_image = numpy.zeros((len(lat_mohid)-1,len(lon_mohid)-1))
        
        my_interpolating_function = RegularGridInterpolator((lat_ds1,lon_ds1),ds1_bnd)
        
        #Interpolate satellite image to MOHID grid
        j = 0
        for lt in lat_mohid:
            i = 0
            for ln in lon_mohid:
                try:
                    interpolated_satellite_image[len(lat_mohid)-j-2][i-1] = my_interpolating_function((lt,ln))
                except:
                    interpolated_satellite_image[len(lat_mohid)-j-2][i-1] = -99
                i = i + 1
            j = j + 1
            
        # register all of the GDAL drivers
        gdal.AllRegister()
        
        # create the output image
        mohid_src = gdal.Open(tif_from_hdf) 
        driver = mohid_src.GetDriver()
        
        values_mohid = mohid_src.GetRasterBand(1).ReadAsArray()
        # Subtract rasters
        difference = interpolated_satellite_image - values_mohid
        difference [difference < -10] = -99
        difference [difference > 50] = -99
        sum_obs = 0
        sum_mod = 0
        n_obs = 0
        n_mod = 0
        for j in range(0, 150):
            for i in range(0, 150):

                if interpolated_satellite_image[i][j] >= 0 and values_mohid[i][j] >= 0:
                
                    sum_obs = sum_obs + interpolated_satellite_image[i][j]
                    n_obs = n_obs + 1
                    
                    sum_mod = sum_mod + values_mohid[i][j]
                    n_mod = n_mod + 1

                        
        avrg_obs = sum_obs/n_obs
        avrg_mod = sum_mod/n_mod               
                        
        num_r2 = 0
        den1 = 0
        den2 = 0
        num_rmse = 0
        num_pbias = 0
        
        for j in range(0, 150):
            for i in range(0, 150):
            
                if interpolated_satellite_image[i][j] >= 0 and values_mohid[i][j] >= 0:
                    num_r2 = num_r2 + (interpolated_satellite_image[i][j] - avrg_obs) * (values_mohid[i][j] - avrg_mod)
                    den1 = den1 + (interpolated_satellite_image[i][j] - avrg_obs)**2
                    den2 = den2 + (values_mohid[i][j] - avrg_mod)**2
                    
                    num_rmse = num_rmse + (interpolated_satellite_image[i][j] - values_mohid[i][j])**2
                    
                    num_pbias = num_pbias + (interpolated_satellite_image[i][j] - values_mohid[i][j])
                    
                    
        den_r2 = (den1**0.5) * (den2**0.5)
        r2 = (num_r2/den_r2)**2
        rmse = (num_rmse/n_mod)**0.5
        pbias = num_pbias/sum_obs
        #print driver
        outDs = driver.Create('differences_' + str(instant) + '.tif', len(lon_mohid), len(lat_mohid), 1, gdal.GDT_Float32)
        if outDs is None:
            print ('Could not create ' + 'differences' + str(instant) + '.tif')
            sys.exit(1)
        
        outBand = outDs.GetRasterBand(1)
        
        # write the data
        outBand.WriteArray(difference, 0, 0)
        
        # flush data to disk, set the NoData value and calculate stats
        outBand.FlushCache()
        outBand.SetNoDataValue(-99)
        
        # georeference the image and set the projection
        outDs.SetGeoTransform(mohid_src.GetGeoTransform())
        outDs.SetProjection(mohid_src.GetProjection())
        #Write the statistical parameters in the file

        #f.write("This is line %d\r\n" % (i+1))
        fi.write(str(instant)+ " ")
        fi.write(str(r2)+ " ")
        fi.write(str(rmse)+ " ")
        fi.write(str(pbias)+ " \n")

        #Close the datasets
        ds1_src = None
        mohid_src = None
        outDs = None
        ds1_bnd = None
        values_mohid = None
        outBand = None
        
        os.remove(ascii_from_hdf)
        os.remove(tif_from_hdf)
        os.remove(satellite_image_wgs)

        f = f + 1