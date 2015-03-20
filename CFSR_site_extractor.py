import netCDF4
import numpy as np
import pandas as pd
import datetime
import sys


def get_T362_lon_lat():
    """
    return the latitude and longitde vectors for the CFSR Guassian
    T362 grid.
    """

    # arbitrarily use soil moisture for July 2008
    URL = Nomads_CFSR_HourlyTS_URL('200807', 'soilm1')
    try:
        nc = netCDF4.Dataset(URL.get_URL())
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        nc.close()
        return(lon, lat)
    except RuntimeError:
        sys.stderr.write('error reading lon, lat')
        sys.stderr.flush()


class SoilSite(object):

    def __init__(self, name, lat, lon):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.data = None  # pandas data frame containing t, soilT, VWC
        self.T362_lonidx = None
        self.T362_latidx = None

    def __str__(self):
        """return a string "Name (lon, lat)" """
        return('{} ({}, {}), ({}, {})'.format(self.name,
                                              self.lon,
                                              self.lat,
                                              self.T362_x,
                                              self.T362_y))

    def get_T362_idx(self, T362_lon, T362_lat):
        """find the x and y indices for the site's latitude and longitude in
        the CFSR Guassian T362 grid.
        """
        self.T362_lonidx = np.searchsorted(T362_lon, self.lon + 360)
        self.T362_latidx = 1 + len(T362_lat) - np.searchsorted(
            T362_lat,
            self.lat,
            sorter=np.arange(len(T362_lat), 0, -1))
        return(self)


class Nomads_CFSR_HourlyTS_URL(object):
    """ build the NOMADS Climate Forecast System Reanalysis hourly
    time series URL for a specified variable and month"""
    def __init__(self, t_str, varname):
        self.t_str = t_str
        self.varname = varname
        self.url_base = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/cfsr1hr/'

    def __str__(self):
        return('{}{}/{}.gdas.{}.grb2'.format(self.url_base,
                                             self.t_str,
                                             self.varname,
                                             self.t_str))

    def get_URL(self):
        return(str(self))


class Nomads_CFSR_HourlyTS(object):
    """reader for NOMADS Climate Forecast System Reanalysis hourly time
    series.  Reads the time series for 5 cm VWC and 5 cm soil T for
    2005 through 2011.

    """
    def __init__(self, year, month, varname_file, varname_nc):
        self.year = year
        self.month = month
        self.varname_file = varname_file
        self.varname_nc = varname_nc

    def read(self, lon_idx, lat_idx):
        """read a grb file from a specified URL
        """
        URL = Nomads_CFSR_HourlyTS_URL('{:04d}{:02d}'.format(self.year,
                                                             self.month),
                                       self.varname_file)
        print(URL)
        try:
            nc = netCDF4.Dataset(URL.get_URL())
        except RuntimeError:
            print('error reading file')

        depth_layer = 0  # these files only have one depth

        t0 = datetime.datetime.now()
        sys.stdout.write('reading {}'.format(self.varname_nc))
        sys.stdout.flush()
        data = nc.variables[self.varname_nc][:, depth_layer, lat_idx, lon_idx]
        sys.stdout.write('...done ({} s)\n'.format(
            (datetime.datetime.now() - t0).seconds))
        sys.stdout.flush()
        self.data = pd.DataFrame(
            data,
            index=pd.DatetimeIndex(freq='1H',
                                   start=datetime.datetime(self.year,
                                                           self.month, 1),
                                   periods=data.size),
            columns=(self.varname_nc,))
        nc.close()
        return(self)


if __name__ == "__main__":

    T362_lon, T362_lat = get_T362_lon_lat()

    names = ['Bondville',
             'Stunt Ranch',
             'Willow Creek',
             'OK',
             'Boyd Deep',
             'Peru']

    lats = [40.0062,  # bondville
            34.1,  # stunt ranch
            45.8060,  # Willow Creek
            36.605,  # OK
            33.648056,  # Boyd Deep
            -12.569167]  # Peru

    lons = [-88.2904,  # bondville
            -118.65,  # stunt ranch
            -90.0798,  # willow creek
            -97.485,  # OK
            -116.376667,  # Boyd Deep
            -70.100111]  # Peru

    sites = [SoilSite(*k) for k in zip(names, lats, lons)]
    sites = [s.get_T362_idx(T362_lon, T362_lat) for s in sites]

    vwc = Nomads_CFSR_HourlyTS(2008, 12,
                               'soilm1',
                               'Volumetric_Soil_Moisture_Content')
    Tsoil = Nomads_CFSR_HourlyTS(2008, 12,
                                 'soilt1',
                                 'Temperature')
    # vwc.read(sites[0].T362_lonidx, sites[0].T362_latidx)
    Tsoil.read(sites[0].T362_lonidx, sites[0].T362_latidx)

    dtidx = pd.DatetimeIndex(freq=pd.tseries.offsets.MonthBegin(),
                             start=datetime.datetime(2008, 1, 1, 0, 0, 0),
                             end=datetime.datetime(2010, 7, 14))

    # time = '200912'
    # url_VWC = Nomads_CFSR_HourlyTS_URL(time, 'soilm1')
    # url_Tsoil = Nomads_CFSR_HourlyTS_URL(time, 'soilt1')

    # ncVWC = read_URL(url_VWC.get_URL())
    # ncTsoil = read_URL(url_Tsoil.get_URL())
