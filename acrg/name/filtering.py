import numpy as np
import numba
import xarray as xr
import acrg.name.name as name

def filtering(fp_data_H,sites,species,filter_types):
    """
    Filters observations in fp_data_H dataset, to account for local influence
    and PBLH thresholds.
    Based on Mark Lunt's code and advice from Alistair Manning.
    Input variables after 'filtering_types' only required for height filtering.
    
    filtering_types = ['localness','pblh','height']
    """

    @numba.jit()
    def assign_localness(fp_data_H_site_fp,wh_rlat,wh_rlon):
        local_sum=np.zeros(fp_data_H_site_fp.shape[2])
        for ti in range(fp_data_H_site_fp.shape[2]):
            fp_data_H_site_local = fp_data_H_site_fp[wh_rlat[0]-2:wh_rlat[0]+3,
                                                     wh_rlon[0]-2:wh_rlon[0]+3,ti]
            local_sum[ti] = np.sum(fp_data_H_site_local)/np.sum(fp_data_H_site_fp[:,:,ti])
        return local_sum

    def define_localness(fp_data_H,site):
        """
        Define the localness of each time point for each site.
        Sum up the 25 grid boxes surrounding the site.
        """
        release_lon=fp_data_H[site].release_lon[0].values
        release_lat=fp_data_H[site].release_lat[0].values
        dlon=fp_data_H[site].lon[1].values-fp_data_H[site].lon[0].values
        dlat=fp_data_H[site].lat[1].values-fp_data_H[site].lat[0].values
        wh_rlon = np.where(abs(fp_data_H[site].lon.values-release_lon) < dlon/2.)[0]
        wh_rlat = np.where(abs(fp_data_H[site].lat.values-release_lat) < dlat/2.)[0]
        fp_data_H_site_fp = fp_data_H[site].fp.values
        local_sum = assign_localness(fp_data_H_site_fp,wh_rlat,wh_rlon)
        local_ds = xr.Dataset({'local_ratio': (['time'], local_sum)},
                                coords = {'time' : (fp_data_H[site].coords['time'])})
        fp_data_H[site] = fp_data_H[site].merge(local_ds)

        return fp_data_H

    def localness_filter(dataset,site,keep_missing=False):
        """
        Subset for times when local influence is below threshold.       
        Local influence expressed as a fraction of the sum of entire footprint domain.
        """
        pc = 0.1      #localness 'ratio' limit
        lr = dataset.local_ratio
        ti = [i for i, local_ratio in enumerate(lr) if local_ratio <= pc]

        if keep_missing is True: 
            mf_data_array = dataset.mf
            dataset_temp = dataset.drop('mf')
            dataarray_temp = mf_data_array[dict(time = ti)]   
            mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            dataset_out = name.combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)] 
            


    n_obs = []
    n_obs_filtered = []

    for i,site in enumerate(sites):
        n_obs.append(fp_data_H[site].mf.values.shape[0])

        if len(filter_types) == len(sites):
            filter_site = filter_types[i]
        else:
            filter_site = filter_types

        if 'local_influence' in filter_site:
            print('Applying localness filtering')
            fp_data_H = define_localness(fp_data_H,site)
            fp_data_H[site] = localness_filter(fp_data_H[site],site)

        if 'BRW wind filt' in filter_site:
            print('BRW wind filt not yet incorporated')

        n_obs_filtered.append(fp_data_H[site].mf.values.shape[0])  

    perc_filtered = np.round((np.array(n_obs) - np.array(n_obs_filtered))/np.array(n_obs)*100,2)
    print(f'% of {species} filtered: {perc_filtered}')   

    return fp_data_H,perc_filtered