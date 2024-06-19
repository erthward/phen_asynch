import matplotlib.pyplot as plt
import geopandas as gpd
import contextily as ctx
import xyzservices
import numpy as np
import os

world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres')) 

def map_taxon(tid, name, fig,
              tiles=False,
              tile_source='OpenTopoMap',
             ):
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)
    fn = f"TID_{tid}_{name}.json"
    obs = gpd.read_file(os.path.join('/media/deth/SLAB/diss/3-phn/inat/obs_data/', fn))
    if not tiles:
        world.plot(ax=ax, zorder=0)
        obs.plot('doy_circ', ax=ax, cmap='twilight', legend=True, zorder=1)
        ax.set_xlim(np.min(obs.geometry.x), np.max(obs.geometry.x))
        ax.set_ylim(np.min(obs.geometry.y), np.max(obs.geometry.y))
    else:
        obs_prj = obs.to_crs(3857)
        obs_prj.plot('doy_circ', ax=ax, cmap='twilight', legend=True, zorder=1)
        ctx.add_basemap(ax=ax, source=xyzservices.providers[tile_source])
        ax.set_xlim(np.min(obs_prj.geometry.x), np.max(obs_prj.geometry.x))
        ax.set_ylim(np.min(obs_prj.geometry.y), np.max(obs_prj.geometry.y))
    ax.set_title(name)


def inspect_inat_phen_maps(res_gdf,
                           tids=None,
                           npeaks=None,
                           show_hist=False,
                           tiles=False,
                          ):
    if npeaks is not None:
        if 'signif_npeaks' not in res_gdf.columns:
            res_gdf.loc[:, 'signif_npeaks'] = [row['npeaks'] if
               row['npeaks_pval']<=0.01 else 0 for i, row in res_gdf.iterrows()]
        tids = res_gdf[res_gdf['signif_npeaks']==npeaks]['tid'].values
    else:
        assert tids is not None
    fig, ax = plt.subplots(1,1)
    for tid in tids:
        name = res_gdf[res_gdf['tid'] == tid]['name'].values[0].replace(' ', '_')
        map_taxon(tid, name, fig, tiles=tiles)
        if show_hist:
            fn = os.path.join('/media/deth/SLAB/diss/3-phn/inat/hist_plots/',
                              f"TID_{tid}_{name}.png")
            cmd = f"feh -ZFd {fn}"
            os.system(cmd)
        input("\n\n<Enter> to proceed to next taxon...\n\n")

