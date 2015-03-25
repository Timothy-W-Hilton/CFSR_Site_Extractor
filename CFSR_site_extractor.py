"""This module reads 5 cm soil temperature (Tsoil) and volumetric
water content (VWC) from the Climate Forecast System (CFS) Version 2
(CFSv2) Reanalysis (CFSR; Saha et al (2010)).  Data are read for the
CFSR grid cells containing specified longitude and latitude
coordinates.  The data are read from the NOAA NOMADS system
(http://nomads.ncdc.noaa.gov/) via web connection.

The data are written to netcdf files, one file per location.

The "__main__" script defines the locations whose grid cells are to be
read, reads the data, writes the netcdf files, and produces a plot of
each time series read.

Author: Timothy W. Hilton, UC Merced  <thilton@ucmerced.edu>
Date: 24 March 2015

REFERENCES

Saha, S., S. Moorthi, H.-L. Pan, X. Wu, J. Wang, S. Nadiga, P. Tripp,
R. Kistler, J. Woollen, D. Behringer, H. Liu, D. Stokes, R. Grumbine,
G. Gayno, J. Wang, Y.-T. Hou, H.-Y. Chuang, H.-M. H. Juang, J. Sela,
M. Iredell, R. Treadon, D. Kleist, P. Van Delst, D. Keyser, J. Derber,
M. Ek, J. Meng, H. Wei, R. Yang, S. Lord, H. Van Den Dool, A. Kumar,
W. Wang, C. Long, M. Chelliah, Y. Xue, B. Huang, J.-K. Schemm,
W. Ebisuzaki, R. Lin, P. Xie, M. Chen, S. Zhou, W. Higgins, C.-Z. Zou,
Q. Liu, Y. Chen, Y. Han, L. Cucurull, R. W. Reynolds, G. Rutledge, and
M. Goldberg (2010), The NCEP Climate Forecast System Reanalysis,
Bulletin of the American Meteorological Society, 91(8), 1015-1057,
doi:10.1175/2010BAMS3001.1.

"""
import netCDF4
import numpy as np
import pandas as pd
import datetime
import sys
import os.path
import matplotlib.pyplot as plt
import glob


def get_T382_land_mask():
    """return the latitude and longitde vectors for the CFSR Guassian T382
    grid (Saha et al (2010), table 1).  The mask is read from the July
    2008 6-houlry forecast for hour 0 (selected arbitrarily).

    From this FAQ (http://rda.ucar.edu/datasets/ds093.0/#docs/FAQs_6hrly.html):

    --------------------------------------------------
    - Where is the land/sea mask in CFSR?

    The land/sea mask is stored in the grid for parameter 'Land
    cover'. Values of '1' indicate land and '0' indicate water. Note:
    There is no land/sea mask data in the hourly time series dataset
    (because the mask does not change with time). You will need to
    obtain the land/sea mask from the main 6-hourly products dataset.
    --------------------------------------------------


    """

    URL = ('http://nomads.ncdc.noaa.gov/thredds/dodsC/modeldata/cmd_flxf'
           '/2008/200807/20080701/flxf01.gdas.2008070100.grb2')
    try:
        nc = netCDF4.Dataset(URL)
        LS_mask = nc.variables['Land_cover_1land_2sea'][...]
        nc.close()
        return(LS_mask)
    except RuntimeError:
        sys.stderr.write('error reading land-sea mask')
        sys.stderr.flush()


def get_T382_lon_lat():
    """return the latitude and longitde vectors for the CFSR Guassian
    T382 grid (Saha et al (2010), table 1).  The coordinates are read
    from the July 2008 soil moisture time series data file (selected
    arbitrarily).
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
    """class to obtain 5 cm CFSR soil moisture and soil temperature data
    from NOMADS (http://nomads.ncdc.noaa.gov/).
    """
    def __init__(self, name, lat, lon):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.data = None  # pandas data frame containing t, soilT, VWC
        self.T382_lonidx = None
        self.T382_latidx = None

    def __str__(self):
        """return a string "Name (lon, lat), (T382_x, T382.y)" T382 is the
        CFSR T382 Gaussian grid.  T382 grid: See Saha et al (2010)
        table 1.

        """
        return('{} ({}, {}), ({}, {})'.format(self.name,
                                              self.lon,
                                              self.lat,
                                              self.T382_lonidx,
                                              self.T382_latidx))

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
        """read hourly 5 cm soil moisture and soil temperature data from
        NOMADS to a pandas data frame with columns Tsoil and VWC and
        indexed by hourly timestamp.  The data frame is placed in the
        object's 'data' field.

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
        """build a netcdf file name for the object's data.  The file is named
        in the format SiteName_Tsoil_VWC.nc, where Site Name is the
        object's name field with spaces removed.

        """

        fname = '{}_Tsoil_VWC.nc'.format(self.name.replace(' ', ''))
        return(fname)

    def data_2_netcdf(self, fname=None):
        """write the object's data field to a netcdf file.  The default file
        name is the result of the get_nc_fname method.  Latitude,
        longitude, T382 grid coordinates, site name, and data citation
        are placed in the netCDF file metadata.

        """
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
        nc.citation = ('Saha, S., S. Moorthi, H.-L. Pan, X. Wu, J. Wang,'
                       ' S. Nadiga, P. Tripp, R. Kistler, J. Woollen, '
                       'D. Behringer, H. Liu, D. Stokes, R. Grumbine, '
                       'G. Gayno, J. Wang, Y.-T. Hou, H.-Y. Chuang, '
                       'H.-M. H. Juang, J. Sela, M. Iredell, R. Treadon, '
                       'D. Kleist, P. Van Delst, D. Keyser, J. Derber, '
                       'M. Ek, J. Meng, H. Wei, R. Yang, S. Lord, '
                       'H. Van Den Dool, A. Kumar, W. Wang, C. Long, '
                       'M. Chelliah, Y. Xue, B. Huang, J.-K. Schemm, '
                       'W. Ebisuzaki, R. Lin, P. Xie, M. Chen, S. Zhou, '
                       'W. Higgins, C.-Z. Zou, Q. Liu, Y. Chen, Y. Han, '
                       'L. Cucurull, R. W. Reynolds, G. Rutledge, and '
                       'M. Goldberg (2010), The NCEP Climate Forecast System '
                       'Reanalysis, Bulletin of the American Meteorological '
                       'Society, 91(8), 1015-1057, '
                       'doi:10.1175/2010BAMS3001.1.')

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
    """build the NOMADS Climate Forecast System Reanalysis hourly time
    series URL for a specified variable and month

    """
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
    series.  Reads the time series for 5 cm VWC and 5 cm soil T for a
    specified name, month, filename, and variable.  For file names see
    http://nomads.ncdc.noaa.gov/docs/CFSR-Timeseries.pdf.

    """
    def __init__(self, year, month, varname_file, varname_nc):
        self.year = year
        self.month = month
        self.varname_file = varname_file
        self.varname_nc = varname_nc

    def read(self, lon_idx, lat_idx):
        """read a single pixel's CFSR time series for the year, month, file,
        and variable specified by the object.

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


class SoilNC(object):
    """make three-panel plot of a site's data, from the netcdf
    file written by SoilSite.data_2_netcdf"""

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


# ======================================================================
def get_data(sites, start_year=2000, end_year=2009):
    """loop through a list of SoilSite objects and read Tsoil and VWC for
    2000 through 2009.

    """
    all_months = pd.DatetimeIndex(freq=pd.tseries.offsets.MonthBegin(),
                                  start=datetime.datetime(start_year, 1, 1),
                                  end=datetime.datetime(end_year, 12, 31))

    for s in sites:
        for this_date in all_months:
            s.get_data(this_date.year, this_date.month)
            s.data_2_netcdf()


def plot_data():
    """Read all the netCDF files generated and plot their data to a time
    series plot, one plot per panel.  This is just a wrapper function
    to make it easier to comment this in or out.

    """
    all_files = glob.glob(os.path.join('/Users/tim/work/Data',
                                       '2015-03-23_CFSR_Soil_Data/*.nc'))
    fig, ax = plt.subplots(nrows=len(all_files), ncols=1, figsize=(8, 24))
    for i, this_file in enumerate(all_files):
        nc = SoilNC(this_file)
        nc.read()
        nc.plot(ax=ax[i])
    fig.savefig('soil_data.pdf')
# ======================================================================


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
    sites = [s.get_T382_idx(T382_lon, T382_lat) for s in sites]

    # Stunt Ranch is 4 miles inland of the Pacific Coast, and happens
    # to fall in a CFSR water cell.  Adjust its longitude index one to
    # the East so that it fall on land.
    sites[1].T382_lonidx += 1

    get_data(sites)
    plot_data()
