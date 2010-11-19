#!/usr/bin/python
def read(name):
    from enthought.tvtk.api import tvtk
    from enthought.mayavi import mlab
    import numpy


    gridreader = tvtk.XMLStructuredGridReader()
    gridreader.file_name = name
    gridreader.update()

    grid = gridreader.output
    data = grid.point_data
    points=grid.points
    phase= numpy.array(data.get_array("Phase"))
    print reversed(grid.dimensions)
    phase_numpy=phase.reshape(grid.dimensions)
    #phase = mlab.pipeline.scalar_field(phase.transpose())
 
    #mlab.pipeline.image_plane_widget(self.vx, plane_orientation='x_axes', slice_index=10)
#mlab.pipeline.image_plane_widget(self.vx, plane_orientation='y_axes', slice_index=10)
#mlab.pipeline.image_plane_widget(self.vx, plane_orientation='z_axes', slice_index=10)
    mlab.contour3d(phase_numpy)
    mlab.show()
    #mlab.axes()

#       mlab.pipeline.vector_cut_plane(self.velocity, mask_points=2, scale_factor=3, plane_orientation='y_axes')
#       mlab.pipeline.vectors(self.velocity, mask_points=20, scale_factor=3.)
    #mlab.outline()


#        # FIXME: Gracefully handle the case when there are no scalar fields.
#        fields = self.sim.output_fields.keys()
#        ffld = fields[0]
#        fields = fields[1:]

#        id.point_data.scalars = self.sim.output_fields[ffld].flatten()
#        id.point_data.scalars.name = ffld

#        for fld in fields:
#            tmp = id.point_data.add_array(self.sim.output_fields[fld].flatten())
#            id.point_data.get_array(tmp).name = fld

#        id.update()

#        for k, v in self.sim.output_vectors.iteritems():
#            if self.sim.grid.dim == 3:
#                tmp = id.point_data.add_array(numpy.c_[v[0].flatten(), v[1].flatten(), v[2].flatten()])
#            else:
#                tmp = id.point_data.add_array(numpy.c_[v[0].flatten(), v[1].flatten(), numpy.zeros_like(v[0].flatten())])
#            id.point_data.get_array(tmp).name = k

#        if self.sim.grid.dim == 3:
#            id.dimensions = list(reversed(self.sim.output_fields[ffld].shape))
#        else:
#            id.dimensions = list(reversed(self.sim.output_fields[ffld].shape)) + [1]
#        w = tvtk.XMLPImageDataWriter(input=id, file_name=('%s%0' + self.digits + 'd.xml') % (self.fname, i))
#        w.write()


if __name__=="__main__":
    read("../phase00000.vts")