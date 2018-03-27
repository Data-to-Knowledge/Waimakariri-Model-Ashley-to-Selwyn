"""
Author: matth
Date Created: 14/03/2017 8:48 AM
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm



class LR(object):
    """
    class to complete a simple linear regression and hold various properties.  easy plot of regression diagnotic plots
    complex numbers are not implemented
    """
    def __init__(self, x, y):
        x = np.array(x)
        y = np.array(y)
        xidx = np.isfinite(x)
        y = y[xidx]
        x = x[xidx]
        yidx = np.isfinite(y)
        y = y[yidx]
        x = x[yidx]
        self.x = x
        self.y = y
        x = sm.add_constant(x)
        self.model = sm.OLS(y,x).fit()

        #below is a bit excessive. it's a legecy issue. you can also call from LR.model
        self.rval = self.model.rsquared
        self.adj_rval = self.model.rsquared_adj
        self.pval = self.model.pvalues
        self.stderr = self.model.mse_resid
        self.corr = np.corrcoef(self.x, self.y)[0, 1]
        self.fitted = self.model.fittedvalues
        self.residuals = self.model.resid
        self.stand_res = self.residuals / self.stderr ** (1 / 2)
        self.slope = self.model.params[1]
        self.intercept = self.model.params[0]
        self.formula = 'y = {}*x + {}'.format(self.slope,self.intercept)

    def predict(self, x):
        if type(x) is int or type(x) is float or type(x) is long or type(x) is np.float64:
            x = [1, x] # add a constant for a single value
            return_only_value = True
        elif len(x) ==1:
            x = [1, x] # add a constant for a single value
        else:
            x = np.array(x)
            x = sm.add_constant(x)
        y_hat = self.model.predict(x)

        if type(x) is int or type(x) is float or type(x) is long or type(x) is np.float64:
            if return_only_value:
                y_hat = y_hat[0]

        return y_hat

    def plot(self, show=True, savepath=None):
        fig,((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,figsize=(18,9))

        # residual vs fitted plot
        ax1.plot(self.fitted, self.residuals, marker='o', linestyle='None')
        ax1.plot(self.fitted, self.residuals * 0, marker='None', linestyle='--', color='r')
        ax1.set_title('Residuals vs. Fitted')
        ax1.set_xlabel('Fitted Values')
        ax1.set_ylabel('Residuals')

        # Normal Q-Q plot
        sm.qqplot(self.residuals, ax=ax2,line='r')
        ax2.set_title('Normal Q-Q')
        ax2.set_xlabel('Theoretical Quantiles')
        ax2.set_ylabel('Standardized Residuals')


        # Scale-Location
        ax3.plot(self.fitted, np.abs(self.stand_res)**(1/2), 'bo')
        ax3.set_title('Scale-Location')
        ax3.set_xlabel('Fitted Values')
        ax3.set_ylabel('|Standardized Residuals|^0.5')
        #could add line from r plot


        # Residuals vs Leverage
        sm.graphics.influence_plot(self.model, ax=ax4)


        # plot data
        ax5.plot(self.x,self.y, 'bo')
        ax5.plot(self.x,self.fitted,'r--')
        ax5.set_title('Fit')
        ax5.set_xlabel('X')
        ax5.set_ylabel('Y')

        ax6.set_xlim(0,5)
        ax6.set_ylim(0,5)
        ax6.text(0.5, 1, 'R2: {}'.format(self.rval))
        ax6.text(0.5, 2, 'adj-R2: {}'.format(self.adj_rval))
        ax6.text(0.5, 3, 'P: int: {}, slope: {}'.format(*self.pval))
        ax6.text(0.5, 3.5, 'params: int: {}, slope: {}'.format(*self.model.params))
        ax6.text(0.5, 4, 'st_err: {}'.format(self.stderr))

        plt.tight_layout()

        if show:
            fig.show()
        if savepath is not None:
            fig.savefig(savepath)
            plt.close()



if __name__ == '__main__':
    x = np.array([1,2,3,4,5,6]) + 10
    y = np.array([1,2,3,4,5,6])
    lm = LR(x,y)
    lm.plot()

    print('test')


