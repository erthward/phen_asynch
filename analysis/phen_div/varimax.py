from numpy import eye, asarray, dot, sum, diag
from numpy.linalg import svd
import matplotlib.pyplot as plt
import numpy as np
from factor_analyzer import FactorAnalyzer


class EOFRotator:
    """
    class to hold raw EOFs and assist in their rotation
    """
    def __init__(self, eofs, x='x', y='y',
                 rotate=True, gamma = 1, q = 20, tol = 1e-6):
        self.eofs = eofs
        self.x = x
        self.y = y
        self.eofs_stack = eofs.stack(spat=[x, y])
        self.n_eofs = self.eofs_stack.shape[0]
        self.n_all_points = self.eofs_stack.shape[1]

        # placeholder for rotated_eofs_indices that are not missing
        self.rotated_eofs_not_missing = None
        # placeholder for stacked eofs subsetted to only non-missing values
        self.eofs_stack_valid = None
        # then fill those placeholders
        self.stack_all_nonmissing_eof_vals()

        # initiate empty array, to be filled with rotated eofs
        self.rotated_eofs = self.eofs * np.nan
        # and same for the eofs_stack
        self.rotated_eofs_stack = self.eofs_stack * np.nan
        # slot to hold the fitted FactorAnalyzer object
        self.fa = None
        # then fill them
        self.rotate=rotate
        self.varimax(rotate=self.rotate)


    def stack_all_nonmissing_eof_vals(self):
        # NOTE: make sure all EOFs are either all valid values or all missing
        # values (no mix-n-match)
        assert (set([*np.unique(np.sum(np.isnan(self.eofs_stack), axis=0))]) ==
                                            set([0, self.n_eofs]))
        # get cols without nans
        not_missing = np.where(np.sum(np.isnan(self.eofs_stack),
                                      axis=0) == 0)[0]

        # store only the non-missing ones, as well as their indices (for later)
        eofs_stack_valid = self.eofs_stack[:, not_missing]
        self.rotated_eofs_not_missing = not_missing
        self.eofs_stack_valid = eofs_stack_valid


    def varimax(self, rotate=True):
        """
        carry out varimax rotation on nxm array stacked EOFs array
        using the factor_analyzer.FactorAnalyzer class
        """
        if rotate:
            rotation = 'varimax'
        else:
            rotation = None
        fa = FactorAnalyzer(n_factors=self.n_eofs,
                            method='principal',
                            rotation=rotation)
        fa.fit(self.eofs_stack_valid[:, :])
        self.rotated_eofs_stack[:,
                                self.rotated_eofs_not_missing] = fa.loadings_.T
        self.rotated_eofs[:,:,:] = self.rotated_eofs_stack.values.reshape(
                                                self.eofs.shape, order='F')
        self.fa = fa


    def plot(self, lyrs=None):
        if lyrs is None:
            lyrs = [*range(self.n_eofs)]
        # get factor variances for labeling
        # TODO: HAD TO HARDCODE ORIGINAL VARIANCES RIGHT NOW BECAUSE CAN'T
        #       FIGURE OUT HOW TO ACCESS RASTER METADATA THROUGH RXR!
        # FIXME!
        orig_var_pcts = [0.77858055, 0.11531084, 0.068770625, 0.037338015]
        orig_var_pcts_norm = [v/np.sum(orig_var_pcts) for v in orig_var_pcts]
        rot_cum_var = self.fa.get_factor_variance()[2]
        rot_var_pct_diffs = rot_cum_var-np.concatenate([[0], rot_cum_var])[:-1]
        rot_var_pcts = rot_var_pct_diffs/rot_cum_var[-1]
        fig = plt.figure(figsize=(8,4*len(lyrs)))
        for i, lyr_n in enumerate(lyrs):
            ax0 = fig.add_subplot(len(lyrs), 2, (i*2+1))
            try:
                self.eofs[lyr_n].plot.imshow(cmap='coolwarm_r', ax=ax0,
                                             add_colorbar=False)
            except AttributeError as e:
                pass
            ax1 = fig.add_subplot(len(lyrs), 2, (i*2+2))
            try:
                self.rotated_eofs[lyr_n].plot.imshow(cmap='coolwarm_r', ax=ax1,
                                                     add_colorbar=False)
            except AttributeError as e:
                pass
            if i == 0:
                ax0.set_title('EOFs')
                if self.rotate:
                    ax1.set_title('varimax-rotated FactorAnalyzer EOFs')
                else:
                    ax1.set_title('unrotated FactorAnalyzer EOFs')
            ax0.set_ylabel('%0.2f%% var' % (100*orig_var_pcts_norm[lyr_n]),
                           fontdict={'fontsize': 18})
            ax1.set_ylabel('%0.2f%% var' % (100*rot_var_pcts[lyr_n]),
                           fontdict={'fontsize': 18})
            ax1.yaxis.set_label_position("right")
            ax1.yaxis.tick_right()
            for ax in [ax0, ax1]:
                ax.set_xticks(())
                ax.set_yticks(())
                ax.set_xlabel('')
                if i>0:
                    ax.set_title('')
        fig.show()
        return fig




# try it out
import rioxarray as rxr
res = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/seasonality/'
                         'results/maps/global_4_EOFs_coswts.tif'))
res_CA = res[:, 400:600, 850:1050]
er = EOFRotator(res_CA)
fig = er.plot()
fig.subplots_adjust(top=0.94, bottom=0.02, left=0.07, right=0.93)
fig.savefig('varimax_rot_small_example_CA.png', dpi=800)
