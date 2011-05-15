#!/usr/bin/python

import numpy as np
from sailfish import lbm, lb_binary
from sailfish import geo

import optparse
from optparse import OptionGroup, OptionParser, OptionValueError

class GeoFE(geo.LBMGeo3D):

    wall_width = 2

    def define_nodes(self):
        hz, hy, hx = np.mgrid[0:self.lat_nz, 0:self.lat_ny, 0:self.lat_nx]

        wall_map = np.logical_or(hy < self.wall_width, hz < self.wall_width)
        self.set_geo(wall_map, self.NODE_WALL)
        self.set_mirror(1, 1)
        self.set_mirror(2, 1)

    def init_fields(self):
        hz, hy, hx = np.mgrid[0:self.lat_nz, 0:self.lat_ny, 0:self.lat_nx]
        width = self.sim.options.iwidth

        self.sim.phi[:] = 1.0
        self.sim.phi[np.logical_and(
                       np.logical_and(hz >= width-1 + self.wall_width, hy >= width-1 + self.wall_width),
                       np.logical_and(hx >= self.lat_nx/3.0 - 1.0, hx<2*self.lat_nx/3.0))] = -1.0
        self.sim.rho[:] = 1.0

class FESim(lb_binary.BinaryFluidFreeEnergy):
    filename = 'fe_separation_3d'

    def __init__(self, geo_class, defaults={}):
        settings = {'verbose': True, 'lat_nx': 750, 'grid': 'D3Q19',
                    'lat_ny': 52, 'lat_nz' : 52,
                    'kappa': 0.04, 'Gamma':1.0, 'A': 0.04,
                    'scr_scale': 1,
                    'tau_a': 2.5, 'tau_b': 0.7, 'tau_phi': 2.0,
                    'bc_wall': 'halfbb', 'bc_wall_grad_order': 1,
                    'periodic_x': True, 'periodic_y': True,'periodic_z': True}
        settings.update(defaults)

        opts = []
        opts.append(optparse.make_option('--force', dest='force', type='float', default=1e-5))
        opts.append(optparse.make_option('--iwidth', dest='iwidth', type='int', default=5))

        lb_binary.BinaryFluidFreeEnergy.__init__(self, geo_class, options=opts, defaults=settings)

        self.add_body_force((self.options.force, 0.0, 0.0), grid=0, accel=False)

        # Use the fluid velocity in the relaxation of the order parameter field,
        # and the molecular velocity in the relaxation of the density field.
        #self.use_force_for_eq(None, 0)
        self.use_force_for_eq(0, 1)


if __name__ == '__main__':
    sim = FESim(GeoFE)
    sim.run()
