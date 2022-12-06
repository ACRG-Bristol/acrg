'''
Make time varying boundary conditions

@author: vf20487

--------------------------------------

Example
-------
# create a BundaryConditions object
bc_obj = BoundaryConditions(filename       = filename,
                            vmr_var        = 'vmr',
                            gph_height_var = 'height',
                            time_coord     = 'time',
                            height_coord   = 'height',
                            species        = 'ch4',
                            domain         = 'EUROPE',
                            start_date     = None)

# create boundary conditions from the input dataset and save to a netcdf file
# all arguments are optional and will default to the ones here
bc_obj.make_bc_file(fp_directory    = None,
                    fp_height_coord = 'height',
                    reverse         = None,
                    convert_units   = True,
                    datasource      = None,
                    out_path        = None,
                    glob_attrs      = {},
                    copy_glob_attrs = False,
                    verbose         = True)

-------

'''
import os
import numpy as np
import xarray as xr
from scipy import interpolate
import matplotlib.pyplot as plt
from acrg.config.paths import Paths
    
class BoundaryConditions:
    def __init__(self, vmr_var, altitude_height_var, dataset=None, filename=None,
                 file_dir=None, time_coord='time', species=None, domain='EUROPE',
                 start_date=None, adjust=None, model=None):
        '''
        vmr_var : str
            The VMR variable name in the nesw dataset
            VMR - vertical mixing ratio
        altitude_height_var : str
            The height variable name in the vmr dataset
            dataset[height_var] will be a 4D array, with time, height,
            lat, and lon coordinates
        dataset (xarray.Dataset, optional)
            dataset containing vmr field from which to create boundary conditions
            if None, filename must be given
        filename (str, optional)
            nc file containing vmr field from which to create boundary conditions
            if None, dataset must be given
        file_dir (str, optional)
            path to directory containing filename
            if None, the filename must include its path
        time_coord (str, optional)
            The time coordinate name, if not 'time'
        species : str
            The gas species of the VMR
        domain : str
            region of study e.g 'EUROPE'
        start_date (str)
            Start date of BC, e.g. '2015-01-01'
        '''
        
        if dataset is None:
            filename = os.path.join(file_dir, filename) if all((file_dir, filename)) else filename
            if filename is None:
                print('Please provide either a dataset or filename containing vmr')
            elif os.path.isfile(filename):
                from acrg.name.name import open_ds
                dataset = open_ds(filename)
            else:
                print(f'Cannot find file: {filename}')
        
        self.data_path    = Paths.data
        self.acrg_path    = Paths.acrg
        
        self.vmr_var      = vmr_var
        self.height_var   = altitude_height_var
        self.time_coord   = time_coord
        self.species      = species
        self.domain       = domain
        self.start_date   = dataset.time.values[0] if start_date is None else start_date
        self.model        = model
        
        # rename latitude and longitude to short version
        if 'latitude'  in dataset: dataset.rename({'latitude'  : 'lat'})
        if 'longitude' in dataset: dataset.rename({'longitude' : 'lon'})

        if adjust is not None:
            print(f'Adjusting vmr by {adjust}')
            dataset[vmr_var] = dataset[vmr_var] + adjust
        self.dataset = dataset
    
    def pressure_to_altitude(self, pressure_var='height', surface_pressure_var='ps', scale_height=8e3,
                             verbose=True):
        '''
        Convert height coordinates in pressure units to height above ground level
        Uses the equation: P = P0 exp(-z / H)
        P = pressure, P0 = surface pressure, z = altitude, H = scale height

        Units of self.dataset[pressure_var] and self.dataset[surface_pressure_var] should match

        Args
            pressure_var (str)
                Variable in self.dataset for the height levels in pressure units
                units of e.g. millibar, bar, hPa, Pa
            surface_pressure_var (str)
                Variable in self.dataset for the pressure at surface level
                units of e.g. millibar, bar, hPa, Pa
            scale_height (int or float)
                The increase in altitude over which the atmospheric pressure
                decreases by a factor of e (approximately 2.718).
                This varies with temperature: H/T = 29.26 m/K
                Units : metres
        '''
        if verbose:
            print(f'Converting height levels from pressure to altitude')
        pressure_ratio = self.dataset[pressure_var] / self.dataset[surface_pressure_var]
        self.dataset['altitude'] = -np.log(pressure_ratio) * scale_height
        self.dataset['altitude'] = self.dataset.altitude.transpose(self.time_coord, pressure_var, 'lat', 'lon')
        self.height_var = 'altitude'

    def make_bc_file(self, fp_directory=None, fp_height_coord='height', reverse=None,
                     convert_units=True, datasource=None, out_path=None, glob_attrs={},
                     copy_glob_attrs=False, verbose=True):
        '''
        Create boundary conditions which match to the NAME grid and save to a netcdf file
        
        Args
            fp_directory (str, optional)
                Path to the NAME footprint files which the input data will be matched to
            fp_height_coord (str, optional)
                Name of the heght coordinate n the footprint file
                Defaults to 'height'
            reverse : (bool/None, optional)
                Whether height values within th dataset input are in reverse order 
                (i.e. nesw["gph"] level 1 values > nesw["gph"] level 2 values).
                Default = None. If this is set to None this will be automatically determined.
            convert_units (bool, optional)
                If True, units will be converted from mol/mol to parts per e.g. ppt, ppb, ppm
                The units to be converted to will be extracted from acrg_species_info.json
            datasource (str, optional)
                source from which data was taken
            out_path : str, optional)
                path to save outputs
            glob_attrs (dict, optional)
                global attributes to add to Dataset with:
                 - key : title of attributes
                 - value : description of attribute
                data will already include
                 - a title with the data source and self.species
                 - the author
                 - the date created
            copy_glob_attrs (bool/list, optional)
                If True, all global attributes from the original dataset will be copied
                to the new file
                Or a list of attributes can be given which will be copied
                If False, none will be copied
            verbose (bool, optional)
                Whether to print any updates
        
        '''
        if verbose: print('\nCutting vmr data to edges\n-------------------------')
        # cut the vmr arry to get just the edges
        self.cut_to_edge(fp_directory = fp_directory, fp_height_coord = fp_height_coord,
                         convert_units = convert_units, verbose = verbose)
        
        if verbose: print('\nInterpolating\n-------------')
        # interpolate to match the NAME height, latitude, and longitude
        self.interpolate_all(reverse=reverse, verbose=verbose)
        
        if verbose: print('\nSaving to netcdf\n----------------')
        # save to a netcdf file with a standardised filename
        self.to_netcdf(datasource=datasource, out_path=out_path, glob_attrs=glob_attrs,
                       copy_glob_attrs=copy_glob_attrs, verbose=verbose)
            
    
    def get_unit_conversion(self, verbose=True):
        '''
        Find the conversion for mol/mol to parts-per required for the gas species
        
        Outputs:
            int
                Adds an integer to the BoundaryConditions object
                Can be accessed using BoundaryConditions_obj.conversion
        '''
        if 'json' not in dir(): import json
        from collections import OrderedDict
        from acrg.convert import concentration
        
        with open(os.path.join(self.acrg_path, "acrg_species_info.json")) as ff:
            species_info=json.load(ff, object_pairs_hook=OrderedDict)
        
        self.conversion = concentration(species_info[self.species.upper()]['units'])
        if verbose: print(f'Unit conversion : {self.conversion}')
    
    def cut_to_edge(self, fp_directory=None, fp_height_coord='height', convert_units=True, verbose=True):
        '''
        Cut the vmr array to the edge of the region needed
        
        Args:
            fp_directory (str, optional)
                Path to the NAME footprint files which the input data will be matched to
                Defualts to <data_path>/LPDM/fp_NAME'
            fp_height_coord (str, optional)
                Name of the heght coordinate n the footprint file
                Defaults to 'height'
            convert_units (bool, optional)
                If True, units will be converted from mol/mol to parts per e.g. ppt, ppb, ppm
                The units to be converted to will be extracted from acrg_species_info.json
            verbose (bool, optional)
                Whether to print any updates
        
        Outputs:
            dict
                Adds a dictionary to the BoundaryConditions object which contans the edge vmrs
                for north, south, east, and west
                Can be accessed using BoundaryConditions_obj.edges
        '''
        from acrg.countrymask import domain_volume
        
        fp_directory = os.path.join(self.data_path, 'LPDM', 'fp_NAME') if fp_directory is None else fp_directory
        self.fp_lat, self.fp_lon, self.fp_height = domain_volume(self.domain, fp_directory=fp_directory)
        
        if convert_units:
            self.get_unit_conversion
            self.dataset[self.species].values *= self.conversion
        
        # Select the gridcells closest to the edges of the  domain and make sure outside of fp
        latlon = {}
        latlon['n'] = (np.abs(self.dataset.coords['lat'].values - max(self.fp_lat))).argmin()
        if self.dataset.coords['lat'].values[latlon['n']] < max(self.fp_lat) and latlon['n'] != 0:
            latlon['n'] -= 1
        
        latlon['s'] = (np.abs(self.dataset.coords['lat'].values - min(self.fp_lat))).argmin()
        if self.dataset.coords['lat'].values[latlon['s']] > min(self.fp_lat) and latlon['s'] != (len(self.dataset.coords['lat'].values)-1):
            latlon['s'] += 1
        
        latlon['e'] = (np.abs(self.dataset.coords['lon'].values - max(self.fp_lon))).argmin()
        if self.dataset.coords['lon'].values[latlon['e']] < max(self.fp_lon) and latlon['e'] != (len(self.dataset.coords['lon'].values)-1):
            latlon['e'] += 1
        
        latlon['w']= (np.abs(self.dataset.coords['lon'].values - min(self.fp_lon))).argmin()
        if self.dataset.coords['lon'].values[latlon['w']] > min(self.fp_lon) and latlon['w'] != 0:
            latlon['w'] -= 1
            
        # Cut to these
        self.edges = {}
        for direction in ['north', 'south']:
            self.edges[direction] = self.dataset.sel(lat = self.dataset.coords['lat'][latlon[direction[0]]],
                                                     lon = slice(self.dataset.coords['lon'][latlon['w']],
                                                                 self.dataset.coords['lon'][latlon['e']])).drop_vars(['lat'])
        for direction in ['east', 'west']:
            self.edges[direction] = self.dataset.sel(lon = self.dataset.coords['lon'][latlon[direction[0]]],
                                                     lat  = slice(self.dataset.coords['lat'][latlon['s']],
                                                                  self.dataset.coords['lat'][latlon['n']])).drop_vars(['lon'])
    
    def interpolate_all(self, reverse=None, verbose=True):
        '''
        Interpolates the data to the NAME heights, latitudes, and longitudes

        Args
            reverse (bool/None, optional)
                Whether height values within th dataset input are in reverse order 
                (i.e. nesw["gph"] level 1 values > nesw["gph"] level 2 values).
                Default = None. If this is set to None this will be automatically determined.
            height_coord (str, optional)
                Used if reverse is not defined to extract appropriate
                height values to compare.
                nesw[height_coord] should be a 1D array of datetimes
                verbose (bool, optional)
                    Whether to print any updates

        Output
            xarray.Dataset
                Adjusts the edges dataset, replacing the 3D VMR arrays for all directions
                with 3D VMR arrays interpolated to the NAME heights, latitudes, and longitudes
                The adjusted dataset can be accessed using BoundaryConditions_obj.edges
        '''
        for dd in ['north', 'south', 'east', 'west']:
            if verbose: print(f'-- Interpolating {dd} boundary vmr --')
            self.interp_height_single(direction=dd, reverse=reverse, verbose=verbose)
            self.interp_latlon_single(direction=dd, verbose=verbose)

    def interp_height_single(self, direction, reverse=None, new_height_coord='height', verbose=True):
        '''
        Interpolates the data to the NAME heights

        Args
            direction (str, optional)
                The compass direction of nesw
                Used for naming the output array : {self.species}_{direction}
            new_height_coord (str, optional)
                Name of the new height variable
                verbose (bool, optional)
            reverse (bool/None, optional)
                Whether height values within is nesw input are in reverse order 
                (i.e. nesw["gph"] level 1 values > nesw["gph"] level 2 values).
                Default = None. If this is set to None this will be automatically determined.
            new_vmr_var (str, optional)
                name for interpolated vmr variable
                defaults to None, inwhich case the original variable name is used
            verbose (bool, optional)
                Whether to print any updates

        Output
            xarray.Dataset
                Adjusts the edges dataset, replacing the 3D VMR array for the given direction
                with a 3D VMR array interpolated to the NAME heights
                Can be accessed using BoundaryConditions_obj.edges
        '''
        # get the axis along which to interpolate
        # either 'longitude' (or 'lon', N or S) or 'latitude' (or 'lat', E or W)
        latorlon = 'lon' if 'lon' in self.edges[direction] else 'lat'
        if verbose: print(f'Interpolating height along {latorlon}')

        # determine whther the heights need to be reversed
        if reverse is None:
            z_check = [float(self.edges[direction][self.height_var][0][hh][0].values) for hh in range(2)]
            reverse = False if z_check[1] >= z_check[0] else True

        # 3D array to fill with interpolated heights
        interp = np.zeros((len(self.edges[direction].coords[self.time_coord]),
                           len(self.fp_height),
                           len(self.edges[direction][latorlon])))

        # loop through the time and latlon coordinates
        # self.edges[direction][height_var] : time, height, latlon
        for i in range(len(self.edges[direction][self.height_var][:,0,0])):
            for j in range(len(self.edges[direction][self.height_var][0,0,:])):
                # get the height and vmr array for all heights
                x = self.edges[direction][self.height_var][i,:,j][::-1] if reverse else \
                    self.edges[direction][self.height_var][i,:,j]
                y = self.edges[direction][self.vmr_var][i,:,j][::-1] if reverse else \
                    self.edges[direction][self.vmr_var][i,:,j]

                # interpolate and create a new array matching NAME heights
                f = interpolate.interp1d(x, y, bounds_error = False, fill_value = np.nan)

                interp[i,:,j] = f(self.fp_height)

        # save interpolated vmr array to new dataset
        self.edges[direction] = xr.Dataset({self.vmr_var : (['time', new_height_coord, latorlon], interp)},
                                           coords = {'time'   : self.edges[direction][self.time_coord].values,
                                                     new_height_coord : self.fp_height,
                                                     latorlon : self.edges[direction][latorlon].values})
        self.edges[direction] = self.edges[direction].ffill(dim=new_height_coord)
        self.edges[direction] = self.edges[direction].bfill(dim=new_height_coord)
    
    def interp_latlon_single(self, direction=None, height_coord='height', reverse=None, verbose=True):
        '''
        Interpolates the data to the NAME latitudes and longitudes

        Args
            direction (str, optional)
                The compass direction of nesw
                Used for naming the output array : {self.species}_{direction}
            height_coord (str, optional)
                Used if reverse is not defined to extract appropriate
                height values to compare.
                nesw[height_coord] should be a 1D array of datetimes
                verbose (bool, optional)
                    Whether to print any updates
            reverse : bool/None
                Whether height values within is nesw input are in reverse order 
                (i.e. nesw["gph"] level 1 values > nesw["gph"] level 2 values).
                Default = None. If this is set to None this will be automatically determined.
                name for interpolated vmr variable
                defaults to None, inwhich case the original variable name is used
            new_vmr_var (str, optional)
                name for interpolated vmr variable
                defaults to None, inwhich case the original variable name is used
            verbose (bool, optional)
                Whether to print any updates

        Output
            xarray.Dataset
                Adjusts the edges dataset, replacing the 3D VMR array for the given direction
                with a 3D VMR array interpolated to the NAME latitudes and longitudes
                Can be accessed using BoundaryConditions_obj.edges
        '''
        # get the axis along which to interpolate
        # either 'longitude' (or 'lon', N or S) or 'latitude' (or 'lat', E or W)
        latorlon = 'lon'       if 'lon'       in self.edges[direction] else \
                   'longitude' if 'longitude' in self.edges[direction] else \
                   'lat'       if 'lat'       in self.edges[direction] else \
                   'latitude'
        if verbose: print(f'Interpolating lat/lon along {latorlon}')

        # determine whther the heights need to be reversed
        if reverse is None:
            ll_check = [self.edges[direction][latorlon].values[ll] for ll in range(2)]
            reverse  = False if ll_check[1] >= ll_check[0] else True

        # 3D array to fill with interpolated heights
        fp_latlon = self.fp_lat if latorlon=='lat' else self.fp_lon
        interp = np.zeros((len(self.edges[direction][self.time_coord]),
                           len(self.fp_height),
                           len(fp_latlon)))

        # loop through height and time
        # self.edges[direction][height_var] : time, height, latlon
        for j in range(len(self.edges[direction][self.vmr_var][0,:,0])):
            for i in range(len(self.edges[direction][self.vmr_var][:,0,0])):
                # get the vmr for all latorlon
                y = self.edges[direction][self.vmr_var][i,j,:][::-1] if reverse else \
                    self.edges[direction][self.vmr_var][i,j,:]
                # get the lat or lon coord
                x = self.edges[direction][latorlon][::-1] if reverse else self.edges[direction][latorlon]

                # interpolate and create a new array matching NAME heights
                f = interpolate.interp1d(x, y, bounds_error = False, fill_value = np.nan)
                interp[i,j,:] = f(fp_latlon)

        # save interpolated vmr array to new dataset
        self.edges[direction] = xr.Dataset({self.vmr_var: (['time', 'height', latorlon], interp)},
                                            coords = {'time'   : self.edges[direction][self.time_coord].values,
                                                      'height' : self.edges[direction][height_coord].values,
                                                      latorlon : fp_latlon})
        self.edges[direction] = self.edges[direction].ffill(dim=latorlon)
        self.edges[direction] = self.edges[direction].bfill(dim=latorlon)
    
    def plot_vmr(self, quantiles=[0.25, 0.5, 0.75], height=500, figsize=(25, 12), fontsize=16, units='ppm',
                 colors=None, cmap='tab10', return_fig=False, show=True, filename=None, savefig_kwargs={}):
        '''
        Plot the VMR at the boundary from the VMR field compared with the boundary condition

        Args
            quantiles (list, optional)
                quantiles along the lat/lon array to plot the VMR timeseries
                list of integers <1
                defaults to quantiles=[0.25, 0.5, 0.75]
            height (int or float, optional)
                height at which to plot the VMR timeseries
                defaults to height=500
            figsize (tuple, optional)
                figure size
                defaults to figsize=(25, 12)
            fontsize (int, optional)
                size of font for the legend & axis labels
                defaults to fontsize=16
            units (str, optional)
                VMR units, used for the axis label
                defaults to units='ppm'
            colors (list or None, optional)
                list of colors to plot the VMR for each quantile
                if None then colors will be selected from cmap
                defaults to colors=None
            cmap (str, optional)
                colormap from which to select colors to plot the VMR for each quantile
                defaults to cmap='tab10'
            return_fig (bool, optional)
                whether to return the fig and axes
            show (bool, optional)
                whether to show the figure
            filename (str, optional)
                filename to save figure
                if not given the figure is not saved
            savefig_kwargs (dict, optional)
                keyword args to use when saving figure

        Output
            matplotlib.Figure and matplotlib axes
                Displays a figure with an axis for each compass direction
        '''
        # get the lat and lon of the domain edges
        bounds = {'north': float(self.fp_lat[-1]), 'east': float(self.fp_lon[-1]),
                  'south': float(self.fp_lat[0]), 'west': float(self.fp_lon[0])}
        # whether to select lat or lon for each boundary
        ll_sel = {'north': 'lon', 'east': 'lat', 'south': 'lon', 'west': 'lat'}
        # use the quantiles to get the lat and lon at which to compare the BCS with the raw field
        latlons = {ll: self.edges[dd][ll].quantile(quantiles).values for ll, dd in {'lon':'north', 'lat':'east'}.items()}
        # create dictionaries for selecting the lat and lon from the bcs and VMR field below
        sel_kwargs_field = {direction: [{'lat':bounds[direction], 'lon':lon, 'method':'nearest'} for lon in latlons['lon']] if ll_sel_dir=='lon' else
                                       [{'lon':bounds[direction], 'lat':lat, 'method':'nearest'} for lat in latlons['lat']]
                            for direction, ll_sel_dir in ll_sel.items()}
        sel_kwargs_bc = {ll: [{ll: latlon, 'height':height, 'method':'nearest'} for latlon in latlon_]
                         for ll, latlon_ in latlons.items()}
        
        # get colors from the colormap unless a list of colors is provided
        colors = [plt.get_cmap(cmap, len(quantiles))(cc) for cc in range(len(quantiles))] if colors is None else colors
        fig, axes = plt.subplots(nrows=4, ncols=1, figsize=figsize, sharex=True, constrained_layout=True)
        for dd, direction in enumerate(bounds.keys()):
            # get the height and VMR arrays for each lat/lon coord
            vmr_quantiles = {quantiles[ss]: self.dataset[self.vmr_var].sel(**sel_kwarg) for ss, sel_kwarg in enumerate(sel_kwargs_field[direction])}
            for qq, quantile in enumerate(quantiles):
                vmr_quantiles[quantile]['height'] = self.dataset[self.height_var].sel(**sel_kwargs_field[direction][qq])[0,:].values

            # plot the boundary condition
            [self.edges[direction][self.vmr_var].sel(**sel_kwarg).plot(ax=axes[dd], color=colors[ss], label=f'{latlons[ll_sel[direction]][ss]:.2f}')
             for ss, sel_kwarg in enumerate(sel_kwargs_bc[ll_sel[direction]])]
            # plot the raw field at the boundary
            [vmr.sel(height=height, method='nearest').plot(ax=axes[dd], ls='--', color=colors[vv])
             for vv, vmr in enumerate(vmr_quantiles.values())]
            # add legend and text to show which boundary this is
            axes[dd].legend(loc='upper right', fontsize=fontsize, title=f'{ll_sel[direction]}:', title_fontsize=fontsize);
            axes[dd].text(0.01, 0.98, f'boundary: {direction}', transform=axes[dd].transAxes, fontsize=fontsize, va='top');
        
        # format the axis labels, title, and ticks
        [ax.set_ylabel(f'VMR, {units}', fontsize=16) for ax in axes];
        [ax.set_xlabel('') for ax in axes];
        [ax.set_title('') for ax in axes];
        [ax.tick_params(labelsize=fontsize) for ax in axes];

        if show:
            fig.show()
        if return_fig:
            return fig, axes
        if filename is not None:
            fig.savefig(filename, **savefig_kwargs)

    def bc_filename(self, from_climatology=False, verbose=True):
        '''
        Create a standardised filename for boundary conditions file
        
        Args
            verbose (bool, optional)
                Whether to print any updates
        
        Returns
            A standardised filename (str)
        '''
        date_str = str(np.datetime64(self.start_date, 'M')).replace('-', '')
        clim_str = "_climatology" if from_climatology else ""
        model_str = f"-{self.model}" if self.model is not None else ""

        self.out_filename = f"{self.species.lower()}{model_str}_{self.domain}_{date_str}{clim_str}.nc"
        if verbose: print(f'Output filename : {self.out_filename}')
    
    def to_netcdf(self, datasource=None, out_path=None, glob_attrs={}, from_climatology=False,
                  copy_glob_attrs=False, verbose=True):
        '''
        Args
            datasource (str, optional)
                source from which data was taken
            out_path (str, optional)
                path to save outputs
            glob_attrs (dict, optional)
                global attributes to add to Dataset with:
                 - key : title of attributes
                 - value : description of attribute
                data will already include
                 - a title with the data source and self.species
                 - the author
                 - the date created
            from_climatology (bool, optional)
                If True, a climatology tag will be added to the name of the output file
            copy_glob_attrs (bool/list, optional)
                If True, all global attributes from the original dataset will be copied
                to the new file
                Or a list of attributes can be given which will be copied
                If False, none will be copied
            verbose (bool, optional)
                Whether to print any updates

        Output
            Creates a netcdf file containing boundary conditions with the heights and latlon
            coordinates interpolated to match the NAME grid coordinates

            Includes DataArrays for the n, s, e, and w boundaries
        '''
        edges = {dd: self.edges[dd].rename({self.vmr_var: f'vmr_{dd[0]}'})
                 for dd in ['north', 'south', 'east', 'west']}
        
        # merge the interpolated n, s, e, & w Datasets
        BC_edges = edges['north'].merge(edges['south']).merge(edges['east']).merge(edges['west'])

        # add global attributes
        BC_edges.attrs['title'] = f'{datasource} {self.species} volume mixing ratios at domain edges'
        BC_edges.attrs['author'] = os.getenv('USER')
        BC_edges.attrs['date_created'] = np.str(np.datetime64('today'))

        # add optional global attributes
        for attr in glob_attrs.keys():
            BC_edges.attrs[attr] = glob_attrs[attr]
        if copy_glob_attrs is True:
            for attr in self.dataset.attrs.keys():
                BC_edges.attrs[attr] = self.dataset.attrs[attr]
        elif isinstance(copy_glob_attrs, list):
            for attr in copy_glob_attrs:
                if attr in self.dataset.attrs.keys():
                    if verbose: print(f'Copying attribute {attr}')
                    BC_edges.attrs[attr] = self.dataset.attrs[attr]
                else:
                    print(f'Could not find {attr} in dataset attributes')

        # create a filename and save dataset to a netcdf file
        out_path = out_path if out_path is not None else \
                   os.path.join(self.data_path, 'LPDM', 'bc', self.domain)
        
        self.bc_filename(from_climatology=from_climatology, verbose=verbose)
        if verbose: print(f'Saving boundary conditions to : {os.path.join(out_path, self.out_filename)}')
        BC_edges.to_netcdf(path = os.path.join(out_path, self.out_filename), mode = 'w')

    