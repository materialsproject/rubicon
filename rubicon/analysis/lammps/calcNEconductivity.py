# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from six.moves import range

__author__ = "mhumbert"


class calcNEconductivity:
    def calcNEconductivity(self, output, molcharge, Lx, Ly, Lz, nummoltype,
                           moltypel, T):
        V = Lx * Ly * Lz * 10 ** -30
        e = 1.60217657e-19
        k = 1.3806488e-23
        NEcond = 0
        for i in range(0, len(moltypel)):
            q = float(molcharge[moltypel[i]])
            if q != 0:
                try:
                    D = float(output['Diffusivity'][moltypel[i]])
                except ValueError:
                    output[
                        'Nernst Einstien Conductivity in S/m'] = 'runtime not long enough'
                    return output
                N = int(nummoltype[i])
                NEcond += N * q ** 2 * D
        NEcond *= e ** 2 / k / T / V

        output['Conductivity']['Nernst_Einstein'] = NEcond

        return output
