from numpy import eye, asarray, dot, sum, diag
from numpy.linalg import svd
import matplotlib.pyplot as plt
import numpy as np
from factor_analyzer import FactorAnalyzer


class EOFRotator:
    """
    class to hold raw EOFs and assist in their rotation
    """
    def __init__(self, eofs, x, y, gamma = 1, q = 20, tol = 1e-6):
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
        self.varimax2()



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


    def varimax(self, gamma = 1, q = 20, tol = 1e-6):
        """
        carry out varimax rotation on nxm array stacked EOFs array,
        given gamma, q, and tol values

        (code stolen from:
            https://stackoverflow.com/questions/44956947/how-to-use-varimax-rotation-in-python)
        """
        # get array shape
        p,k = self.eofs_stack_valid.shape
        # get starting rotatin matrix as a kxk identity matrix
        R = eye(k)
        # start d at 0
        d=0
        # for q times...
        for i in range(q):
            # store previous d value
            d_old = d
            # dot product of stacked EOFs and R
            # NOTE: this just returns the eofs on first loop,
            #       but not in subsequent loops
            Lambda = dot(self.eofs_stack_valid, R)
            # ???
            u,s,vh = (svd(dot(self.eofs_stack_valid.T,asarray(Lambda)**3 - (gamma/p) *
                              dot(Lambda, diag(diag(dot(Lambda.T,Lambda)))))))
            # get new rotation matrix
            R = dot(u,vh)
            d = sum(s)
            # break if the ratio between previous and current 
            if d/d_old < tol: break
        # store the output in the rotated eofs stack
        # (in the columns corresponding to the non-missing values)
        self.rotated_eofs_stack[:, self.rotated_eofs_not_missing] = dot(
                                                      self.eofs_stack_valid, R)
        # and reshape back to original
        self.rotated_eofs[:,:,:] = self.rotated_eofs_stack.values.reshape(
                                                 self.eofs.shape, order='F')


    def varimax2(self):
        """
        carry out varimax rotation on nxm array stacked EOFs array
        using the factor_analyzer.FactorAnalyzer class
        """
        fa = FactorAnalyzer(n_factors=self.n_eofs,
                            method='principal',
                            rotation='varimax')
        fa.fit(self.eofs_stack_valid[:, :])
        self.rotated_eofs_stack[:,
                                self.rotated_eofs_not_missing] = fa.loadings_.T
        self.rotated_eofs[:,:,:] = self.rotated_eofs_stack.values.reshape(
                                                self.eofs.shape, order='F')
        self.fa = fa


    def plot(self, lyr_n=0):
        fig = plt.figure(figsize=(8,4))
        ax0 = fig.add_subplot(121)
        self.eofs[lyr_n].plot.imshow(cmap='coolwarm_r', ax=ax0)
        ax0.set_title('EOFs')
        ax1 = fig.add_subplot(122)
        self.rotated_eofs[lyr_n].plot.imshow(cmap='coolwarm_r', ax=ax1)
        ax1.set_title('varimax-rotated EOFs')
        fig.show()
        return fig
