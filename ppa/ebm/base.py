from abc import ABCMeta, abstractmethod
from ppa.ebm.pathway import AbsPathway

import pandas as pd
from scipy.integrate import solve_ivp
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['AbsModel']


class AbsModel(metaclass=ABCMeta):
    def __init__(self, t0, t1, dt, pp: AbsPathway):
        self.Year0 = t0
        self.Year1 = t1
        self.dYear = dt
        self.YearRange = [t0, t1]
        self.YearSeries = np.linspace(t0, t1, int((t1 - t0) / dt) + 1)
        self.PatientPathway = pp

    @abstractmethod
    def get_y0(self, t, pars) -> np.ndarray:
        pass

    @abstractmethod
    def __call__(self, t, y, pars) -> np.ndarray:
        pass

    @abstractmethod
    def measure(self, t, y, pars):
        pass

    def simulate(self, pars, y0=None):
        if 'pp' not in pars and self.PatientPathway is not None:
            pars['pp'] = self.PatientPathway(pars)

        if y0 is None:
            y0 = self.get_y0(self.Year0, pars)
        ys = solve_ivp(self, self.YearRange, y0=y0, args=(pars, ), dense_output=True)
        ms = pd.DataFrame([self.measure(t, ys.sol(t), pars) for t in self.YearSeries]).set_index('Time')
        return ys, ms


if __name__ == '__main__':
    import matplotlib.pylab as plt

    class SIR(AbsModel):
        def get_y0(self, t, pars) -> np.ndarray:
            return np.array([990, 10, 0])

        def __call__(self, t, y, pars) -> np.ndarray:
            n = y.sum()
            s, i, _ = y
            foi = pars['beta'] * i / n
            return np.array([
                - foi * s,
                foi * s - pars['gamma'] * i,
                pars['gamma'] * i
            ])

        def measure(self, t, y, pars):
            n = y.sum()
            s, i, _ = y
            foi = pars['beta'] * i / n

            return {
                'Time': t,
                'I': i / n,
                'Inc': foi * s / n
            }

    sir = SIR(0, 30, dt=0.5, pp=None)
    _, ms = sir.simulate({'beta': 1.5, 'gamma': 0.2})
    ms.plot()
    plt.show()
