from paraview.simple import *
import time

reader=XMLStructuredGridReader(FileName="/home/shurik/Documents/Projects/binary_microchannel_3d/Temp/phase040000.vts")
reader.UpdatePipeline()

contour=Contour()
contour.ContourBy="Phase"
contour.Isosurfaces=[0.0]
Show(contour)

SetActiveSource(reader)
plane = Slice(SliceType="Plane" )
plane.SliceType.Origin = [1, 25, 1000]
plane.SliceType = "Plane"
plane.SliceOffsetValues=[0.0]


surf=SurfaceVectors()
surf.SelectInputVectors=['Velocity']

seedpoints = MaskPoints(Input=surf)
seedpoints.MaximumNumberofPoints = 200
seedpoints.Random=1

stream=StreamTracerWithCustomSource(Source=seedpoints,Input=surf)

stream.Vectors = ['POINTS', 'Velocity']
stream.MaximumStreamlineLength = 1499.0
#stream.UpdatePipeline()

representation=Show(stream)
PVLookupTable = GetLookupTableForArray( "Velocity", 3, RGBPoints=[-0.00036603630720772327, 0.23000000000000001, 0.29899999999999999, 0.754, 0.023360469851033357, 0.70599999999999996, 0.016, 0.14999999999999999], VectorMode='Component', ColorSpace='Diverging', ScalarRangeInitialized=1.0, VectorComponent=2 )

Velocity_PiecewiseFunction = CreatePiecewiseFunction()

representation.ColorArrayName = 'Velocity'
representation.LookupTable = PVLookupTable
representation.ColorAttributeType = 'POINT_DATA'

Render()

view=GetRenderView()


#view.CameraFocalPoint = [380.15643757605432, -69.418119533743777, 349.24287984226908]
#view.CameraClippingRange = [468.00899956856784, 1917.5261259682597]
numpoints=41
startz=0
endz=1245
totalmovie=10.0

for counter,z in enumerate(range(startz,endz+1,(endz-startz)/(numpoints-1))): 
    view.CameraPosition = [300, 0, z]
    #view.ViewTime =20*totalmovie/(numpoints-1)*counter
    view.CameraFocalPoint = [0, 0, z]
    view.CameraParallelScale = 1.0
    time.sleep(0.2)  
    Render()


#view = GetRenderView()
#view.CameraPosition = [3282.1772141912752, 1495.7211400491408, -294.03426471994089]
#view.ViewTime = 0.1111111111111111
#view.CameraFocalPoint = [1.0, 25.4999, 749.49800000000005]
#view.CameraParallelScale = 1.0
#Render()

#view.ViewTime = 0.22222222222222221
#view.CameraPosition = [1883.4852227731355, 917.32527989081495, -2356.6886617055666]
#Render()

#view.ViewTime = 0.33333333333333331
#view.CameraPosition = [-391.03280116543021, -76.073594268930776, -2969.9686349067579]
#Render()

#view.ViewTime = 0.44444444444444442
#view.CameraPosition = [-2488.063895595943, -1024.5125548640813, -1846.0831563954853]
#Render()

#view.ViewTime = 0.55555555555555558
#view.CameraPosition = [-3419.7655558312335, -1481.4875578502229, 499.78060845066011]
#Render()

#view.ViewTime = 0.66666666666666663
#view.CameraPosition = [-2745.4130067558831, -1230.8347695012826, 2957.9539895293642]
#Render()

#view.ViewTime = 0.77777777777777779
#view.CameraPosition = [-792.16616591387765, -395.01090938796142, 4380.1394105686531]
#Render()

#view.ViewTime = 0.88888888888888884
#view.CameraPosition = [1532.9160579051013, 637.75484483781963, 4110.3029967989933]
#Render()

#view.ViewTime = 1.0
#view.CameraPosition = [3146.9335318226335, 1386.6477729024177, 2264.1134566579599]

#Render()

