#!/usr/bin/python

import numpy as np
from sailfish import lbm
from sailfish import geo

import optparse
from optparse import OptionGroup, OptionParser, OptionValueError

class GeoFE(geo.LBMGeo2D):

    wall_width = 2

    def define_nodes(self):
        hy, hx = np.mgrid[0:self.lat_ny, 0:self.lat_nx]
        wall_map = np.logical_or(hy < self.wall_width, hy >= self.lat_ny-self.wall_width)
        self.set_geo(wall_map, self.NODE_WALL)

    def init_fields(self):
        hy, hx = np.mgrid[0:self.lat_ny, 0:self.lat_nx]
        width = self.sim.options.iwidth

        self.sim.phi[:] = 1.0
        self.sim.phi[np.logical_and(
                       np.logical_and(hx >= (self.lat_nx-1)/3.0,
                                      hx <= 2*(self.lat_nx-1)/3.0),
                       np.logical_and(hy >= width-1 + self.wall_width, hy <= self.lat_ny-width-self.wall_width))] = -1.0
        self.sim.rho[:] = 1.0

class FESim(lbm.BinaryFluidFreeEnergy):
    filename = 'fe_seperation_2d'

    def __init__(self, geo_class, defaults={}):
        settings = {'verbose': True, 'lat_nx': 3001,
                    'lat_ny': 200+2, 'grid': 'D2Q9',
                    'kappa': 0.04, 'Gamma': 1.0, 'A': 0.04,
                    'scr_scale': 1,
                    'tau_a': 2.5, 'tau_b': 0.7, 'tau_phi': 1.0,
                    'bc_wall': 'halfbb',
                    'periodic_x': True, 'periodic_y': True}
        settings.update(defaults)

        opts = []
        opts.append(optparse.make_option('--force', dest='force', type='float', default=0.0))
        opts.append(optparse.make_option('--iwidth', dest='iwidth', type='int', default=5))

        lbm.BinaryFluidFreeEnergy.__init__(self, geo_class, options=opts, defaults=settings)

        self.add_body_force((self.options.force, 0.0), grid=0, accel=False)

        # Use the fluid velocity in the relaxation of the order parameter field,
        # and the molecular velocity in the relaxation of the density field.
        self.use_force_for_eq(None, 0)
        self.use_force_for_eq(0, 1)


if __name__ == '__main__':
    sim = FESim(GeoFE)
    sim.run()
