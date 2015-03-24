"""REFERENCES


"""
import netCDF4
import numpy as np
import pandas as pd
import datetime
import sys
import os.path
import matplotlib.pyplot as plt
import glob


def get_T382_lon_lat():
    """
    return the latitude and longitde vectors for the CFSR Guassian
    T382 grid (Saha et al (2010), table 1).
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
        self.T382_lonidx = None
        self.T382_latidx = None

    def __str__(self):
        """return a string "Name (lon, lat), (T382_x, T382.y)" T382 is the CFSR T382 Gaussian grid"""
        return('{} ({}, {}), ({}, {})'.format(self.name,
                                              self.lon,
                                              self.lat,
                                              self.T382_x,
                                              self.T382_y))

    def get_T382_idx(self, T382_lon, T382_lat):
        """find the x and y indices for the site's latitude and longitude in
        the CFSR Guassian T382 grid.
        """
        self.T382_lonidx = np.searchsorted(T382_lon, self.lon + 360)
        self.T382_latidx = 1 + len(T382_lat) - np.searchsorted(
            T382_lat,
            self.lat,
            sorter=np.arange(len(T382_lat), 0, -1))
        return(self)

    def get_data(self, year, month):
        """
        """
        vwc = Nomads_CFSR_HourlyTS(year, month,
                                   'soilm1',
                                   'Volumetric_Soil_Moisture_Content')
        Tsoil = Nomads_CFSR_HourlyTS(year, month,
                                     'soilt1',
                                     'Temperature')
        vwc.read(self.T382_lonidx, self.T382_latidx)
        Tsoil.read(self.T382_lonidx, self.T382_latidx)

        df = pd.merge(vwc.data, Tsoil.data,
                      left_index=True, right_index=True)

        if self.data is None:
            self.data = df
        else:
            self.data = pd.concat([self.data, df])
        return(self)

    def get_nc_fname(self):
        """build a netcdf file name for the object's data"""
        fname = '{}_Tsoil_VWC.nc'.format(self.name.replace(' ', ''))
        return(fname)

    def data_2_netcdf(self, fname=None):

        NC_DOUBLE = 'd'  # netCDF4 specifier for NC_DOUBLE datatype
        NC_INT64 = 'i8'  # netCDF4 specifier for NC_DOUBLE datatype

        if fname is None:
            fname = self.get_nc_fname()
        nc = netCDF4.Dataset(fname, 'w')
        nc.createDimension('time', self.data.shape[0])
        nc.createVariable(varname='tstamp',
                          datatype=NC_INT64,
                          dimensions=('time'))
        nc.createVariable(varname='VWC',
                          datatype=NC_DOUBLE,
                          dimensions=('time'))
        nc.createVariable(varname='Tsoil',
                          datatype=NC_DOUBLE,
                          dimensions=('time'))

        nc.site = self.name
        nc.site_latitude = self.lat
        nc.site_longitude = self.lon
        nc.T382_x = self.T382_lonidx
        nc.T382_y = self.T382_latidx
        nc.description = ('NCEP Climate Forecast System Reanalysis hourly 5'
                          ' cm soil temperature (Tsoil) and volumetric water'
                          ' content (VWC).  T382_x and T382_y are the site''s'
                          ' x and y indices in the CFSR T382 grid (see '
                          'Saha, et al (2010), table 1).')
        nc.citation = 'Saha, S., S. Moorthi, H.-L. Pan, X. Wu, J. Wang, S. Nadiga, P. Tripp, R. Kistler, J. Woollen, D. Behringer, H. Liu, D. Stokes, R. Grumbine, G. Gayno, J. Wang, Y.-T. Hou, H.-Y. Chuang, H.-M. H. Juang, J. Sela, M. Iredell, R. Treadon, D. Kleist, P. Van Delst, D. Keyser, J. Derber, M. Ek, J. Meng, H. Wei, R. Yang, S. Lord, H. Van Den Dool, A. Kumar, W. Wang, C. Long, M. Chelliah, Y. Xue, B. Huang, J.-K. Schemm, W. Ebisuzaki, R. Lin, P. Xie, M. Chen, S. Zhou, W. Higgins, C.-Z. Zou, Q. Liu, Y. Chen, Y. Han, L. Cucurull, R. W. Reynolds, G. Rutledge, and M. Goldberg (2010), The NCEP Climate Forecast System Reanalysis, Bulletin of the American Meteorological Society, 91(8), 1015-1057, doi:10.1175/2010BAMS3001.1.'

        nc.variables['tstamp'].description = 'time stamp'
        nc.variables['tstamp'].units = 'hours since {}'.format(
            self.data.index[0])
        nc.variables['tstamp'][:] = (
            (self.data.index - self.data.index[0])
            .astype('timedelta64[h]').values
        )

        nc.variables['VWC'].description = '5 cm soil volumetric water content'
        nc.variables['VWC'].units = 'fraction'
        nc.variables['VWC'][:] = self.data[
            'Volumetric_Soil_Moisture_Content'].values

        nc.variables['Tsoil'].description = '5 cm soil temperature'
        nc.variables['Tsoil'].units = 'Kelvins'
        nc.variables['Tsoil'][:] = self.data[
            'Temperature'].values
        nc.close()


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
        try:
            data = nc.variables[self.varname_nc][:, depth_layer,
                                                 lat_idx, lon_idx]
            sys.stdout.write('...done ({} s, {})\n'.format(
                (datetime.datetime.now() - t0).seconds,
                datetime.datetime.now()))
            sys.stdout.flush()
        except UnboundLocalError:
            print "unable to read {}".format(URL)
            self.data = pd.DataFrame(None)
            return(self)
        self.data = pd.DataFrame(
            data,
            index=pd.DatetimeIndex(freq='1H',
                                   start=datetime.datetime(self.year,
                                                           self.month, 1),
                                   periods=data.size),
            columns=(self.varname_nc,))
        nc.close()
        return(self)


class soil_nc(object):
    """make three panel plot of a site's data, from the netcdf
    file written by data_2_netcdf"""

    def __init__(self, fpath):
        self.fpath = fpath

    def read(self):
        nc = netCDF4.Dataset(self.fpath, 'r')
        t = nc.variables['tstamp'][:]
        t_dt64 = (map(lambda this_t: np.timedelta64(this_t, 'h'), t) +
                  np.datetime64('2000-01-01T00:00:00Z'))
        VWC = nc.variables['VWC'][:]
        Tsoil = nc.variables['Tsoil'][:]

        self.data = pd.DataFrame(index=t_dt64,
                                 data=np.dstack((VWC, Tsoil)).squeeze(),
                                 columns=('VWC', 'Tsoil'))
        self.site = nc.site
        self.lat = nc.site_latitude
        self.lon = nc.site_longitude
        nc.close()
        return(self)

    def plot(self, ax=None):
        # plt.style.use('fivethirtyeight')
        df = self.data
        df.Tsoil -= 273.15
        ax = self.data.plot(secondary_y=['VWC'], colormap='Dark2', ax=ax)

        ax.set_xlabel('date')
        ax.set_ylabel('T$_{soil}$ (K)')
        ax.right_ax.set_ylabel('VWC')
        ax.set_title('{} ({} E, {} N)'.format(self.site, self.lon, self.lat))
        return(ax)


if __name__ == "__main__":

    T382_lon, T382_lat = get_T382_lon_lat()

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
    # sites = [s.get_T382_idx(T382_lon, T382_lat) for s in sites]

    # all_months = pd.DatetimeIndex(freq=pd.tseries.offsets.MonthBegin(),
    #                               start=datetime.datetime(2000, 1, 1),
    #                               end=datetime.datetime(2009, 12, 31))

    # for s in sites:
    #     for this_date in all_months:
    #         s.get_data(this_date.year, this_date.month)
    #         s.data_2_netcdf()

    all_files = glob.glob(os.path.join('/Users/tim/work/Data',
                                       '2015-03-23_CFSR_Soil_Data/*.nc'))
    fig, ax = plt.subplots(nrows=len(all_files), ncols=1, figsize=(8, 24))
    for i, this_file in enumerate(all_files):
        nc = soil_nc(this_file)
        nc.read()
        nc.plot(ax=ax[i])
    fig.savefig('soil_data.pdf')
