


# Read in soil temperature data (assumes this is always there)
# ref: http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
soilData = self.header[3]
self.nSoil = int(soilData[1])           # Number of ground temperature depths
self.Tsoil = utilities.zeros(self.nSoil,12)  # nSoil x 12 matrix for soil temperture (K)
self.depth = utilities.zeros(self.nSoil,1)   # nSoil x 1 matrix for soil depth (m)

# Read monthly data for each layer of soil from EPW file
for i in xrange(self.nSoil):
    self.depth[i][0] = float(soilData[2 + (i*16)]) # get soil depth for each nSoil
    # Monthly data
    for j in xrange(12):
        self.Tsoil[i][j] = float(soilData[6 + (i*16) + j]) + 273.15 # 12 months of soil T for specific depth




# Define Road (Assume 0.5m of asphalt)
kRoad = ipd['kRoad']                 # road pavement conductivity (W/m K)
cRoad = ipd['cRoad']                 # road volumetric heat capacity (J/m^3 K)

emis = 0.93
asphalt = Material(kRoad,cRoad,'asphalt')
road_T_init = 293.
road_horizontal = 1
road_veg_coverage = min(vegCover/(1-bldDensity),1.) # fraction of surface vegetation coverage

# define road layers
road_layer_num = int(math.ceil(d_road/0.05))
thickness_vector = map(lambda r: 0.05, range(road_layer_num)) # 0.5/0.05 ~ 10 x 1 matrix of 0.05 thickness
material_vector = map(lambda n: asphalt, range(road_layer_num))

road = Element(alb_road,emis,thickness_vector,material_vector,road_veg_coverage,road_T_init,road_horizontal)


print road
print '--'
print 'soil properties'
print 'soilnum', self.nSoil
print 'soilDepth: ',
pprint.pprint(self.depth)
print 'soilT summer: ',
pprint.pprint( map(lambda t: round(sum(t[2:6])/float(len(t[2:6]))-273.15,2), self.Tsoil))
print '\n'

#testing alt
#self.depth[0][0] = 0.4
#min_depth_diff = float('Inf')
#for i in xrange(self.nSoil):
#    curr_depth_diff = abs(sum(road.layerThickness) - self.depth[i][0])
#    if min_depth_diff >= curr_depth_diff:
#        min_depth_diff = curr_depth_diff
#        soilindex1 = i
