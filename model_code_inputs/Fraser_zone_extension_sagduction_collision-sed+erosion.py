# %% [markdown]
# # Lithostatic model with extension, melting and crust formation
#
# The module can be imported as follows:

# %%
from underworld import UWGeodynamics as GEO
from underworld import visualisation as vis


import underworld as uw

import scipy
from scipy.interpolate import interp1d

import numpy as np
import os
import math



from underworld import function as fn

# %%
render = False

# %%
GEO.rcParams['swarm.particles.per.cell.2D'] = 30

GEO.rcParams['popcontrol.particles.per.cell.2D'] = 30
GEO.rcParams['popcontrol.aggressive']= False
### half the PPC ? 
GEO.rcParams['popcontrol.max.splits'] = 15

# GEO.rcParams

# %% [markdown]
# ### Scaling

# %%
u = GEO.UnitRegistry

half_rate = 1.0 * u.centimeter / u.year
model_length = 660e3 * u.meter
surfaceTemp = 293.15 * u.degK
baseLABTemp = 1673.15 * u.degK
bodyforce = 3300 * u.kilogram / u.metre**3 * 9.81 * u.meter / u.second**2

KL = model_length
Kt = KL / half_rate
KM = bodyforce * KL**2 * Kt**2
KT = (baseLABTemp - surfaceTemp)

GEO.scaling_coefficients["[length]"] = KL
GEO.scaling_coefficients["[time]"] = Kt
GEO.scaling_coefficients["[mass]"]= KM
GEO.scaling_coefficients["[temperature]"] = KT

# %%
GEO.rcParams['temperature.SIunits'] = u.celsius
GEO.rcParams['pressureField.SIunits'] = u.kilobar

# %%
H       = 130. * u.kilometer ### depth of box
L       = 660  * u.kilometer

H_air   = 10 * u.kilometer
H_uCrust = 20 * u.kilometer    ### depth of upper crust
H_lCrust = 40 * u.kilometer    ### depth of lower crust
H_LAB   = 120 * u.kilometer    ### depth of LAB
T0      = 293.15 * u.degK      ### surface temp

rho_new_crust = 3000. * u.kilogram / u.metre**3

T_Moho = (700 + 273.15) * u.degK
T_LAB  = (1375 + 273.15) * u.degK ### LAB temp

basal_HF = 25e-3*u.watt/u.meter**2

adiabatic_gradient = 0.5*u.kelvin/u.kilometer
Tz = T_LAB + (adiabatic_gradient) * (H - H_LAB).to(u.kilometer) ### Bottom of box temp
Tz

D_surf = 150 *u.meter**2/u.year

# %% [markdown]
# #### Define the external geometry
#

# %%
extension_rate = 1*u.centimeter/u.year
extension_duration = 17*u.megayear
total_extension = (extension_rate * extension_duration).to(u.kilometer)
total_extension

# %%
static_duration = 40 * u.megayear

collision_duration = 10 * u.megayear

# %%
Model = GEO.Model(elementRes=(330, 70 ), 
                  minCoord=(0. * u.kilometer, -H), 
                  maxCoord=(L, H_air), 
                  gravity=(0.0, -9.81 * u.meter / u.second**2))

# %%
# Model.outputDir = (f"Fraser_zone_model-"
#                    f"D_LAB={int(H_LAB.m)}km-"
#                    f"T_LAB={int(T_LAB.m - 273.15)}C"
#                    f"-t_ext={int(extension_duration.m)}Myr-ext_r={int(extension_rate.m)}cmyr-1"
#                    f"-t_static={int(static_duration.m)}Myr-ext"
#                    f"-t_col={int(collision_duration.m)}Myr-col_r={int(extension_rate.m)}cmyr-1"
#                    f"-D_surf={int(D_surf.m)}m2yr-1"
#                    f"-nprocs={uw.mpi.size}")

Model.outputDir = (f"Fraser_zone_model-"
                   f"rho_new_crust={int(rho_new_crust.m)}kgm-3-"
                   f"D_ucrust={int(H_uCrust.m)}km-"
                   f"D_Moho={int(H_lCrust.m)}km-"
                   f"T_Moho={int(T_Moho.m - 273.15)}C-"
                   f"D_LAB={int(H_LAB.m)}km-"
                   f"T_LAB={int(T_LAB.m - 273.15)}C-"
                   f"v_ext={extension_rate.m}cmyr-1-"
                   f"t_ext={extension_duration.m}Myr-"
                   # f"t_ext={int(extension_duration.m)}Myr-ext_r={int(extension_rate.m)}cmyr-1-"
                   # f"t_static={int(static_duration.m)}Myr-ext-"
                   # f"t_col={int(collision_duration.m)}Myr-col_r={int(extension_rate.m)}cmyr-1-"
                   f"D_surf={int(D_surf.m)}m2yr-1-"
                   f"xres={Model.elementRes[0]}-"
                   f"yres={Model.elementRes[1]}-"
                   f"nprocs={uw.mpi.size}")

Model.outputDir

# %% [markdown]
# #### Add some Materials
#

# %%
Model.diffusivity = 1e-6 * u.metre**2 / u.second 
Model.capacity    = 1000. * u.joule / (u.kelvin * u.kilogram)

# %%
air = Model.add_material(name="Air", shape=GEO.shapes.Layer(top=Model.top, bottom=0.0 * u.kilometer))
# stickyAir = Model.add_material(name="StickyAir", shape=GEO.shapes.Layer(top=air.bottom, bottom= 0.0 * u.kilometer))

sediment = Model.add_material(name="sediment")
crust1 = Model.add_material(name="Crust1")
crust2 = Model.add_material(name="Crust2")

crust3 = Model.add_material(name="Crust3")
crust4 = Model.add_material(name="Crust4")

mantleLithosphere = Model.add_material(name="MantleLithosphere", shape=GEO.shapes.Layer(top=-H_lCrust, bottom=-H_LAB))
mantle = Model.add_material(name="Mantle", shape=GEO.shapes.Layer(top=mantleLithosphere.bottom, bottom=Model.bottom))




newCrust        = Model.add_material(name="newCrust") ## crust formed from mantle melt

depletedMantle  = Model.add_material(name="newCrust") ## mantle melt has come from

partiallyMoltenMantle = Model.add_material(name="moltenMantle") ## mantle containing melt





# %%
### Make a layered crust 5 km thick
sin_function = np.sign(np.sin(GEO.dimensionalise(Model.swarm.data[:,1], u.kilometer)/(1.6 * u.kilometer)))

Model.materialField.data[(sin_function>0)  &  (Model.swarm.data[:,1] < GEO.nd(0*u.kilometer)) & (Model.swarm.data[:,1] >= GEO.nd(-H_lCrust))] = crust3.index
Model.materialField.data[(sin_function<0)  &  (Model.swarm.data[:,1] < GEO.nd(0*u.kilometer)) & (Model.swarm.data[:,1] >= GEO.nd(-H_lCrust))] = crust4.index


Model.materialField.data[(sin_function>0)  &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer)) & (Model.swarm.data[:,1] >= GEO.nd(-H_uCrust))] = crust1.index
Model.materialField.data[(sin_function<0)  &  (Model.swarm.data[:,1] < GEO.nd(0.*u.kilometer)) & (Model.swarm.data[:,1] >= GEO.nd(-H_uCrust))] = crust2.index


# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200, 400))
    Fig.Points(Model.swarm, Model.materialField, fn_size=4.0)
    Fig.show()

# %% [markdown]
# #### Add some material properties

# %%
air.diffusivity = 1.0e-6 * u.metre**2 / u.second
#. stickyAir.diffusivity = 1.0e-6 * u.metre**2 / u.second

air.capacity = 100. * u.joule / (u.kelvin * u.kilogram)
#. stickyAir.capacity = 100. * u.joule / (u.kelvin * u.kilogram)

# %%
air.density = 1. * u.kilogram / u.metre**3
#. stickyAir.density = 1. * u.kilogram / u.metre**3

sediment.density = GEO.LinearDensity(2600. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

crust1.density = GEO.LinearDensity(2750. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
crust2.density = GEO.LinearDensity(2750. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

crust3.density = GEO.LinearDensity(2950. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

crust4.density = GEO.LinearDensity(2950. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)

newCrust.density = GEO.LinearDensity( rho_new_crust , thermalExpansivity=3e-5 / u.kelvin)


mantleLithosphere.density = GEO.LinearDensity(3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
mantle.density = GEO.LinearDensity(3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)


depletedMantle.density = GEO.LinearDensity(3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)
partiallyMoltenMantle.density = GEO.LinearDensity(3300. * u.kilogram / u.metre**3, thermalExpansivity=3e-5 / u.kelvin)


# %%
crust1.radiogenicHeatProd = 1.5 * u.microwatt / u.meter**3
crust2.radiogenicHeatProd = 1.5 * u.microwatt / u.meter**3

crust3.radiogenicHeatProd = 0.5 * u.microwatt / u.meter**3
crust4.radiogenicHeatProd = 0.5 * u.microwatt / u.meter**3

sediment.radiogenicHeatProd = 1.5 * u.microwatt / u.meter**3

mantleLithosphere.radiogenicHeatProd = 0.02 * u.microwatt / u.meter**3

### need to check values
newCrust.radiogenicHeatProd = 0.2 * u.microwatt / u.meter**3 ### need to check
depletedMantle.radiogenicHeatProd = 0.001 * u.microwatt / u.meter**3 ### need to check

partiallyMoltenMantle.radiogenicHeatProd = 0.005 * u.microwatt / u.meter**3 ### need to check

# %% [markdown]
# ### Define Viscosities


# %%
### Strong crust from Ranelli 1995
Diabase_Dislocation_Ranalli_1995 = GEO.ViscousCreep(preExponentialFactor=2.0e-4/u.megapascal**3.4/u.second,
                                                      stressExponent=3.4,
                                                      activationVolume=0.,
                                                      activationEnergy=260 * u.kilojoules/u.mole,
                                                      waterFugacity=0.0,
                                                      grainSize=0.0,
                                                      meltFraction=0.,
                                                      grainSizeExponent=0.,
                                                      waterFugacityExponent=0.,
                                                      meltFractionFactor=0.0,
                                                      f=1.0)

### Quartzite from Ranalli 1995
Quartzite_Dislocation_Ranalli_1995 = GEO.ViscousCreep(preExponentialFactor=6.7e-6/u.megapascal**2.4/u.second,
                                                      stressExponent=2.4,
                                                      activationVolume=0.,
                                                      activationEnergy=156 * u.kilojoules/u.mole,
                                                      waterFugacity=0.0,
                                                      grainSize=0.0,
                                                      meltFraction=0.,
                                                      grainSizeExponent=0.,
                                                      waterFugacityExponent=0.,
                                                      meltFractionFactor=0.0,
                                                      f=1.0)


# %%
rh = GEO.ViscousCreepRegistry()

# %%
Model.minViscosity = 1e19 * u.pascal * u.second
Model.maxViscosity = 1e24 * u.pascal * u.second

air.viscosity                = 1e19 * u.pascal * u.second
#. stickyAir.viscosity          = 1e20 * u.pascal * u.second

sediment.viscosity            = 0.5*rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
crust1.viscosity              = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
crust2.viscosity              = rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

crust3.viscosity              = rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998
crust4.viscosity              = rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998

mantleLithosphere.viscosity   = rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993
mantle.viscosity              = rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993

# %%
newCrust.viscosity              = rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998
depletedMantle.viscosity        = rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993
partiallyMoltenMantle.viscosity = rh.Dry_Olivine_Dislocation_Karato_and_Wu_1993
# partiallyMoltenMantle.viscosity = 1e18 * u.pascal * u.second

# %% [markdown]
# ### Define Plasticity

# %%
pl = GEO.PlasticityRegistry()

# %%
# crustPlasticity = GEO.DruckerPrager(
#                                      cohesion=10. * u.megapascal,
#                                      cohesionAfterSoftening=1. * u.megapascal,
#                                      frictionCoefficient = 0.4,
#                                      frictionAfterSoftening = 0.04,
#                                      epsilon1=0.0,
#                                      epsilon2=0.5)

# mantlePlasticity = GEO.DruckerPrager(
#                                      cohesion=20. * u.megapascal,
#                                      cohesionAfterSoftening=10 * u.megapascal,
#                                      frictionCoefficient = 0.6,
#                                      frictionAfterSoftening = 0.06,
#                                      epsilon1=0.0,
#                                      epsilon2=0.5)

# %%
sediment.plasticity           = pl.Rey_et_al_2014_UpperCrust

crust1.plasticity             = pl.Rey_et_al_2014_UpperCrust
crust2.plasticity             = pl.Rey_et_al_2014_UpperCrust

crust3.plasticity             = pl.Rey_et_al_2014_UpperCrust
crust4.plasticity             = pl.Rey_et_al_2014_UpperCrust

mantleLithosphere.plasticity  = pl.Rey_et_al_2014_LithosphericMantle
mantle.plasticity             = pl.Rey_et_al_2014_LithosphericMantle

# %%
newCrust.plasticity               = pl.Rey_et_al_2014_UpperCrust
depletedMantle.plasticity         = pl.Rey_et_al_2014_LithosphericMantle
partiallyMoltenMantle.plasticity  = pl.Rey_et_al_2014_LithosphericMantle

# %%
# pl.Rey_et_al_2014_LithosphericMantle

# %% [markdown]
# ## Temperature Boundary Conditions

# %%
Model.set_temperatureBCs(top=T0, 
                        #bottom=Tz,
                         # bottom=-basal_HF,
                         materials=[(air, T0)]) #, (stickyAir, T0)])

# Model.set_heatFlowBCs(bottom=(-basal_HF,
#                              mantle))

# %% [markdown]
# ## Velocity Boundary Conditions

# %%
extension_condition = [(Model.y < GEO.nd(0. * u.kilometre), GEO.nd(extension_rate)),
                       (True, GEO.nd(extension_rate) + Model.y * (GEO.nd((-2. * extension_rate) / GEO.nd(Model.maxCoord[1]))))]


right_wall_vbc = fn.branching.conditional(extension_condition)

# %%
# side_outflux = extension_rate * (Model.maxCoord[1] - Model.minCoord[1])
# side_outflux.to_base_units()



# %%
# top_influx = side_outflux / Model.maxCoord[0]
# top_influx.to_base_units()

# %%
Model.set_velocityBCs(left=[0., None],
                       right=[right_wall_vbc, None],
                       top=[None, 0.], )
                       # bottom=GEO.LecodeIsostasy(reference_mat=mantle, average=False))

# %% [markdown]
# ## Initialize plastic strain

# %%
import numpy as np

def gaussian(xx, centre, width):
    return ( np.exp( -(xx - centre)**2 / width**2 ))

# Set the seed to give the same sequence of random numbers
np.random.seed(0)


Model.plasticStrain.data[:] = 0.
maxDamage = 0.7
Model.plasticStrain.data[:] = maxDamage * np.random.rand(*Model.plasticStrain.data.shape[:])
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,0], GEO.nd(180*u.kilometer), GEO.nd(100.0 * u.kilometer))
Model.plasticStrain.data[:,0] *= gaussian(Model.swarm.particleCoordinates.data[:,1], GEO.nd(-H_lCrust) , GEO.nd(15.0 * u.kilometer))

Model.plasticStrain.data[:,0][Model.swarm.data[:,1] > 0.] = 0.



# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200, 400))
    Fig.Points(Model.swarm, Model.plasticStrain, fn_size=4.0)
    Fig.show()

# %% [markdown]
# ## Passive tracers

# %% [markdown]
# ### Interface tracers

# %%
import numpy as np

npoints = 1000 # This is the number of points used to define the surface
coords = np.ndarray((npoints, 2))

coords[:, 0] = np.linspace(GEO.nd(Model.minCoord[0]), GEO.nd(Model.maxCoord[0]), npoints)
coords[:, 1] = 0.

Model.add_passive_tracers(name="Surface", vertices=coords, advect=False)

coords[:, 1] -= GEO.nd(H_lCrust)
Model.add_passive_tracers(name="Moho", vertices=coords, advect=False)

# %%
# Model.Surface_tracers.allow_parallel_nn = True
# Model.Moho_tracers.allow_parallel_nn = True

# %% [markdown]
# ### Grid Tracers

# %%
pts = GEO.circles_grid(radius=2.0*u.kilometer, 
                    minCoord=[Model.minCoord[0], -H_lCrust], 
                    maxCoord=[Model.maxCoord[0], 0.*u.kilometer])

Model.add_passive_tracers(name="FSE_Crust", vertices=pts)

pts = GEO.circles_grid(radius=2.0*u.kilometer, 
                    minCoord=[Model.minCoord[0], mantleLithosphere.bottom], 
                    maxCoord=[Model.maxCoord[0], mantleLithosphere.top])

Model.add_passive_tracers(name="FSE_lithosphere", vertices=pts)

# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200,400), title="Material Field", quality=2)
    # Fig.Points(Model.Central_tracers, pointSize=10, colour="blue")
    Fig.Points(Model.Surface_tracers, pointSize=2.0)
    Fig.Points(Model.Moho_tracers, pointSize=2.0, colour="red")
    Fig.Points(Model.FSE_Crust_tracers, pointSize=2.0)
    Fig.Points(Model.FSE_lithosphere_tracers, pointSize=2.0)
    Fig.Points(Model.swarm, Model.materialField, fn_size=2.0)
    Fig.show()

# %% [markdown]
# ### Add in melting

# %%
''' melting can occur in the: mantle, lithosphere, newCrust and partiallyMoltenMantle'''
solidii = GEO.SolidusRegistry()
# crust_solidus = solidii.Crustal_Solidus

liquidii = GEO.LiquidusRegistry()
# crust_liquidus = liquidii.Crustal_Liquidus

mantle_solidus = solidii.Mantle_Solidus
mantle_liquidus = liquidii.Mantle_Liquidus

mantleLithosphere.add_melt_modifier(mantle_solidus, mantle_liquidus, 
                         latentHeatFusion=400.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.5,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.2,
                         viscosityChangeX2 = 0.3,
                         viscosityChange = 1e-3
                        )  

mantle.add_melt_modifier(mantle_solidus, mantle_liquidus, 
                         latentHeatFusion=400.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.5,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.2,
                         viscosityChangeX2 = 0.3,
                         viscosityChange = 1e-3
                        ) 

# newCrust.add_melt_modifier(mantle_solidus, mantle_liquidus, 
#                          latentHeatFusion=250.0 * u.kilojoules / u.kilogram / u.kelvin,
#                          meltFraction=0.,
#                          meltFractionLimit=0.3,
#                          meltExpansion=0.13, 
#                          viscosityChangeX1 = 0.2,
#                          viscosityChangeX2 = 0.3,
#                          viscosityChange = 1e-3
#                         ) 

partiallyMoltenMantle.add_melt_modifier(mantle_solidus, mantle_liquidus, 
                         latentHeatFusion=400.0 * u.kilojoules / u.kilogram / u.kelvin,
                         meltFraction=0.,
                         meltFractionLimit=0.5,
                         meltExpansion=0.13, 
                         viscosityChangeX1 = 0.2,
                         viscosityChangeX2 = 0.3,
                         viscosityChange = 1e-3
                        ) 

# %%
# import pyMelt as m
# def calculate_mantle_melt_using_pyMelt():
#     ### using pyMelt
#     lz = m.lithologies.katz.lherzolite()
    
#     temp_swarm = GEO.dim(Model.temperature.evaluate(Model.swarm.data), u.degC).m
#     pressure_swarm = GEO.dim(Model.pressureField.evaluate(Model.swarm.data), u.pascal).m / 1e9
    
#     ### only find the swarm markers where melt occurs, rather than looping through all particles
#     ind = np.argwhere( (Model.materialField.data[:,0] == mantleLithosphere.index) |  (Model.materialField.data[:,0] == mantle.index) | (Model.materialField.data[:,0] == partiallyMoltenMantle.index) )
    
#     meltFraction = np.zeros_like(temp_swarm)
    
#     for i in ind:
#         meltFraction[i] = lz.F(temp_swarm[i], pressure_swarm[i])
    
#     Model.meltField.data[:,] = meltFraction

# %%
# # ### Add in cell centroids using a swarm
# centroids = uw.swarm.Swarm(Model.mesh)
# layout = uw.swarm.layouts.PerCellGaussLayout(centroids, gaussPointCount=1)
# # centroids.populate_using_layout(layout)

# ### Add in grid using coordinates, useful if the mesh is deformed
# # centroids = uw.swarm.Swarm(Model.mesh)
# centroids.add_particles_with_coordinates(grid_coords)


# centroids.allow_parallel_nn = True

# centroid_meltFrac = centroids.add_variable('double', 1)

# %%
# fig = vis.Figure(figsize=(800,400))
# fig.append( vis.objects.Points(swarm=centroids, pointSize=1, colourBar=False) )
# fig.append( vis.objects.Mesh(Model.mesh))
# fig.show()

# %% [markdown]
# # Compute initial condition

# %%

m0 = (T_Moho - T0) / H_lCrust
m1 = (T_LAB - T_Moho) / (H_LAB - H_lCrust )
m2 = (Tz - T_LAB) / (H - H_LAB)

crust_geotherm        = GEO.nd(m0) * (-Model.y) + GEO.nd(T0)
lithosphere_geotherm  = GEO.nd(m1) * (-Model.y - GEO.nd(H_lCrust)) + GEO.nd(T_Moho)
mantle_geotherm       = GEO.nd(m2) * (-Model.y - GEO.nd(H_LAB)) + GEO.nd(T_LAB)



conditions = [(Model.y <= GEO.nd(-H_LAB), mantle_geotherm),
              (Model.y <= GEO.nd(-H_lCrust), lithosphere_geotherm),
              (Model.y <= 0., crust_geotherm),
              (True, GEO.nd(T0))]

geotherm_fn = fn.branching.conditional(conditions)




# %%
### calculate our own lithostatic pressure
Model.init_model(temperature=geotherm_fn, pressure='lithostatic')

# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200,400), title="Pressure Field (MPa)", quality=3)
    Fig.Surface(Model.mesh, GEO.dimensionalise(Model.pressureField, u.megapascal))
    Fig.show()

# %% [markdown]
# ### Calculate lithostatic pressure

# %%
# elementType = "Q1"
# lp_mesh = uw.mesh.FeMesh_Cartesian(elementType=elementType,
#                              elementRes=Model.elementRes,
#                              minCoord=(GEO.nd(Model.minCoord[0]), GEO.nd(Model.minCoord[1])),
#                              maxCoord=(GEO.nd(Model.maxCoord[0]), GEO.nd(Model.maxCoord[1])),
#                              periodic=Model.periodic)

# lithoPressureField        = lp_mesh.add_variable( nodeDofCount=1 )
# velocityField             = lp_mesh.add_variable( nodeDofCount=2 )
# velocityField.data[:]= 0.

# # Boundary conditions
# topWall    = p_mesh.specialSets["MaxJ_VertexSet"]
# bottomWall = p_mesh.specialSets["MinJ_VertexSet"]

# lithoPressureBC = uw.conditions.DirichletCondition( variable        = lithoPressureField, 
#                                                indexSetsPerDof = ( topWall) )




# %%
# #### Add in cell centroids using a swarm for the LP solver
# centroids = uw.swarm.Swarm(Model.mesh)
# layout = uw.swarm.layouts.PerCellGaussLayout(centroids, gaussPointCount=1)
# centroids.populate_using_layout(layout)

# lithoPressureSolver = uw.systems.SteadyStateDarcyFlow(velocityField=velocityField,
#                                             pressureField=lithoPressureField, 
#                                             fn_diffusivity = 1.,
#                                             conditions=[lithoPressureBC],
#                                             fn_bodyforce=Model._buoyancyFn, 
#                                             voronoi_swarm=centroids)



# %%
# LPsolver = uw.systems.Solver(lithoPressureSolver)
# LPsolver.solve()
# if GEO.nProcs == 1 & render:
#     Fig = vis.Figure(figsize=(1200,400), title="Pressure Field (MPa)", quality=3)
#     Fig.Surface(Model.mesh, GEO.dimensionalise(lithoPressureField, u.megapascal))
#     Fig.show()

# %%
# ### project onto UWGeo pressure field
# Model.pressureField.data[...] = lithoPressureField.evaluate(Model.mesh.subMesh)

# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200,400), title="Pressure Field (MPa)", quality=3)
    Fig.Surface(Model.mesh, GEO.dimensionalise(Model.pressureField, u.megapascal))
    Fig.show()

# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200,400), title="Temperature Field (Celsius)", quality=3)
    Fig.Surface(Model.mesh, GEO.dimensionalise(Model.temperature, u.kelvin), colours="coolwarm")
    Fig.show()

# %%
# if GEO.nProcs == 1 & render:
#     import matplotlib.pyplot as plt
#     y_prof = np.arange(10, -100, -0.1)
#     x_prof = np.zeros_like(y_prof)+(Model.maxCoord[0]/2).m
#     profile = np.column_stack([x_prof,y_prof])
#     GEO.nd(profile * u.kilometer)
#     plt.plot(Model.pressureField.evaluate(GEO.nd(profile * u.kilometer)), y_prof)
#     plt.plot(Model.lithostatic_pressureField.evaluate(GEO.nd(profile * u.kilometer)), y_prof)
    
#     plt.plot(lithoPressureField.evaluate(GEO.nd(profile * u.kilometer)), y_prof, ls=':')

# %% [markdown]
# ## Basic Analysis of the initial set-up
#
# Let's analyze the pressure and temperature field by creating a vertical profile at the center of the model.

# %%
# Only run this when in serial. Will fail in parallel

if GEO.nProcs == 1 & render:

    moho_average_temperature = Model.temperature.evaluate(Model.Moho_tracers).mean()
    moho_average_temperature = GEO.dimensionalise(moho_average_temperature, u.degC)

    print("Average Temperature at Moho: {0:5.0f}".format(moho_average_temperature))

    distances, temperature = GEO.extract_profile(Model.temperature, line = [(180.* u.kilometer, 0.), (180.* u.kilometer, Model.bottom)])
    distances, pressure = GEO.extract_profile(Model.pressureField, line = [(180.* u.kilometer, 0.), (180.* u.kilometer, Model.bottom)])

    Fig, (ax1, ax2) = plt.subplots(1,2,figsize=(15,7))
    ax1.plot(GEO.dimensionalise(temperature, u.degK), GEO.dimensionalise(distances, u.kilometer))
    ax1.set_xlabel("Temperature in Kelvin")
    ax1.set_ylabel("Depth in kms")
    ax1.set_ylim(120, 0)
    ax1.set_title("Temperature profile")

    ax2.plot(GEO.dimensionalise(pressure, u.megapascal), GEO.dimensionalise(distances, u.kilometer))
    ax2.set_xlabel("Pressure in megapascal")
    ax2.set_ylabel("Depth in kms")
    ax2.set_title("Pressure profile")
    ax2.set_ylim(120, 0)
    plt.show()

# %% [markdown]
# In a similar fashion, one can extract a vertical profile of the viscosity field.

# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200,400), title="Viscosity Field (Pa.s)", quality=3)
    Fig.Points(Model.swarm, 
               GEO.dimensionalise(Model.viscosityField, u.pascal * u.second),
               logScale=True,
               fn_size=3.0)
    Fig.show()

# %%
# Quartzite_Dislocation_Ranalli_1995
# Diabase_Dislocation_Ranalli_1995

# crust1.viscosity              = Quartzite_Dislocation_Ranalli_1995#rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995
# crust2.viscosity              = Quartzite_Dislocation_Ranalli_1995#rh.Wet_Quartz_Dislocation_Gleason_and_Tullis_1995

# crust3.viscosity              = Diabase_Dislocation_Ranalli_1995#rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998
# crust4.viscosity              = Diabase_Dislocation_Ranalli_1995#rh.Dry_Maryland_Diabase_Dislocation_Mackwell_et_al_1998

# %%
# Only run this when in serial. Will fail in parallel

if GEO.nProcs == 1 & render:

    import matplotlib.pyplot as plt

    distances, viscosities = GEO.extract_profile(Model.projViscosityField, 
                                                 line = [(Model.maxCoord[0] / 2, 0.), (Model.maxCoord[0] / 2, Model.bottom)])
    distances, stresses = GEO.extract_profile(Model.projStressField, 
                                             line = [(Model.maxCoord[0] / 2, 0.), (Model.maxCoord[0] / 2, Model.bottom)])


    Fig, (ax1) = plt.subplots(1,1,figsize=(7,7))
    ax1.plot(GEO.dimensionalise(viscosities, u.pascal * u.second), GEO.dimensionalise(distances, u.kilometer), label='rheology')
    # ax1.plot(GEO.dimensionalise(viscosities1, u.pascal * u.second), GEO.dimensionalise(distances, u.kilometer), label='ranalli')

    ax1.legend()
    ax1.set_xlabel("Viscosity in Pa.s")
    ax1.set_ylabel("Depth in kms")
    ax1.set_ylim(100,-10)
    ax1.set_title("Viscosity profile")

# %% [markdown]
# And the melt field

# %%
if GEO.nProcs == 1 & render:

    T_s = mantle_solidus.temperature(pressure)
    T_l = mantle_liquidus.temperature(pressure)
    depths = distances

    import pylab as plt

    fig = plt.figure(figsize=(8,4))
    plt.plot(GEO.dimensionalise(T_s, u.degC),GEO.dimensionalise(depths, u.kilometer), label="T_s")
    plt.plot(GEO.dimensionalise(T_l, u.degC),GEO.dimensionalise(depths, u.kilometer), label="T_l")
    plt.plot(GEO.dimensionalise(temperature, u.degC),GEO.dimensionalise(depths, u.kilometer), label="Geotherm")
    plt.xlabel("Temperature (C)")
    plt.ylabel("Depth (km)")
    plt.legend()
    plt.ylim(110, 0)
    plt.show()


# %%
def gather_data(val, masked_value=-9999999, bcast=False):

    """
    gather values on root (bcast=False) or all (bcast = True) processors
    Parameters:
        vals : Values to combine into a single array on the root or all processors
        masked_value : value to use as mask value
        bcast : broadcast array/value to all processors

    returns:
        val_global : combination of values form all processors

    """

    comm = uw.mpi.comm
    rank = uw.mpi.rank
    size = uw.mpi.size

    dtype = val.dtype


    ### make sure all data comes in the same order
    # with uw.mpi.call_pattern(pattern="sequential"):
    if len(val > 0):
        val_local = np.ascontiguousarray(val.copy(), dtype=dtype)
    else:
        val_local = np.ascontiguousarray([masked_value], dtype=dtype)


    comm.barrier()

    ### Collect local array sizes using the high-level mpi4py gather
    sendcounts = np.array(comm.gather(len(val_local), root=0))

    if rank == 0:
        val_global = np.zeros((sum(sendcounts)), dtype=dtype)

    else:
        val_global = None

    comm.barrier()

    ## gather x values, can't do them together
    comm.Gatherv(sendbuf=val_local, recvbuf=(val_global, sendcounts), root=0)

    comm.barrier()

    if uw.mpi.rank == 0:
        ### remove rows with NaN
        val_global = val_global[val_global != masked_value]  

    comm.barrier()

    if bcast == True:
        #### make available on all processors
        val_global = comm.bcast(val_global, root=0)

    comm.barrier()

    return val_global


# %%
def vertical_advection(tracers):
    '''
    Function to advect passive tracers in the vertical direction only
    '''
    tracer_vel_x = gather_data( Model.velocityField.evaluate(tracers.data)[:,0] )
    tracer_vel_y = gather_data( Model.velocityField.evaluate(tracers.data)[:,1] )
    tracer_coords_x = gather_data( tracers.data[:,0] )
    tracer_coords_y = gather_data( tracers.data[:,1] )

    tracer_f = None 
    
    if uw.mpi.rank == 0:
        ### advect tracers
        tracer_coords = np.column_stack([tracer_coords_x, tracer_coords_y])
        tracer_vel = np.column_stack([tracer_vel_x, tracer_vel_y])
        new_coords = tracer_coords[:,] + (Model.dt.value*tracer_vel[:,])
        tracer_f = interp1d(new_coords[:,0], new_coords[:,1], fill_value="extrapolate", kind='linear')

    tracer_f = uw.mpi.comm.bcast(tracer_f, root=0)

    uw.mpi.comm.Barrier()

    with tracers.deform_swarm():
        tracers.particleCoordinates.data[:,1] = tracer_f(tracers.particleCoordinates.data[:,0])
        
def advect_tracers_vertically():
    vertical_advection(Model.Surface_tracers)
    vertical_advection(Model.Moho_tracers)


# %%
def diffusive_surface():
    ### function to diffuse the surface to represent erosion and sedimentation

    comm = uw.mpi.comm

    
    surface_tracers_x = gather_data( Model.Surface_tracers.data[:,0] )
    surface_tracers_y = gather_data( Model.Surface_tracers.data[:,1] )

    f1 = None

    if uw.mpi.rank == 0:
        ### sort gathered coords
        surface_coords = np.column_stack([surface_tracers_x, surface_tracers_y])
        surface_coords = surface_coords[np.argsort(surface_coords[:,0])]
        # surface_coords = surface_coords[surface_coords[:, 0].argsort()]

        dx = np.diff(surface_coords[:,0]).min()

        nd_D = GEO.nd(D_surf.to_base_units())

        x_nd = surface_coords[:,0]
        z_nd = surface_coords[:,1]

        ### diffuse the surface
        diffusion_CFL_dt = 0.2 * (dx**2/nd_D)

        timestep = Model.dt.value

        nts = math.ceil(timestep/diffusion_CFL_dt)
            
        surf_dt = (timestep / nts)

        print('SP total time:', GEO.dim(timestep, u.year), 'timestep:', GEO.dim(surf_dt, u.year), 'No. of its:', nts, flush=True)


        ### Basic Hillslope diffusion
        for i in range(nts):
            qs = -nd_D * np.diff(z_nd)/np.diff(x_nd)
            dzdt = -np.diff(qs)/np.diff(x_nd[:-1])


            z_nd[1:-1] += dzdt*surf_dt
    

        ### creates function for the new surface that has eroded, to be broadcast back to nodes
        f1 = interp1d(x_nd, z_nd, fill_value='extrapolate', kind='linear')

    comm.barrier()

    '''broadcast the new surface'''
    ### broadcast function for the surface
    f1 = comm.bcast(f1, root=0)

    with Model.Surface_tracers.deform_swarm():
        Model.Surface_tracers.particleCoordinates.data[:,1] = f1(Model.Surface_tracers.particleCoordinates.data[:,0])

    comm.barrier()

    ### update the time of the sediment and air material as sed & erosion occurs
    if Model.timeField:
        ### Set newly deposited sediment time to 0 (to record deposition time)
        Model.timeField.data[ (Model.swarm.data[:,1] < f1(Model.swarm.data[:,0])) & (Model.materialField.data[:,0] == air.index) ] = 0.
        ### reset air material time back to the model time
        Model.timeField.data[ (Model.swarm.data[:,1] > f1(Model.swarm.data[:,0])) & (Model.materialField.data[:,0] != air.index) ] = Model.timeField.data.max()

    '''Erode surface/deposit sed based on the surface'''
    ### update the material on each node according to the spline function for the surface
    Model.materialField.data[(Model.swarm.data[:,1] > f1(Model.swarm.data[:,0])) & (Model.materialField.data[:,0] != air.index) ] = air.index
    Model.materialField.data[(Model.swarm.data[:,1] < f1(Model.swarm.data[:,0])) & (Model.materialField.data[:,0] == air.index) ] = sediment.index


    
    
Model.post_solve_functions['surface_processes'] = diffusive_surface


# %%
### get centroid of mesh cells
centroid = uw.swarm.Swarm(Model.mesh)
layout = uw.swarm.layouts.PerCellGaussLayout(centroid, gaussPointCount=1)
centroid.populate_using_layout(layout)


# %%
### Amount of curent melt
Model.meltField 
### Amount of total melt produced
Model.cumulativeMelt = Model.add_swarm_variable('cumulativeMelt', projected='mesh', restart_variable=True) ### used to track the total amount of melt on a particle

### amount of melt actually extracted at current timestep
Model.meltExtract = Model.add_swarm_variable('meltExtract', projected='mesh', restart_variable=True) ### used to track the current melting event


# %%
def update_melt_fractions(meltField, cumulativeMelt, min_melt_in_rock, min_extractable_melt_fraction, max_extractable_melt_fraction):
    # (2) Compute how much melt will be extracted and crustal thickness
    total_melt = Model.meltField.data[:,0]
    extracted_melt = Model.cumulativeMelt.data[:,0]
    # (2a) Compute the amount of melt that can be extracted
    available_melt = total_melt - extracted_melt
    available_melt[available_melt < 0] = 0
    
    # (2b) Extract melt if it exceeds a threshold
    melt_above_threshold = available_melt > min_extractable_melt_fraction
    
    melt_for_crust_generation = np.zeros_like(available_melt)
    melt_for_crust_generation[melt_above_threshold] = available_melt[melt_above_threshold] - min_melt_in_rock
    
    # (2c) Do not extract melt from regions that have already reached the maximum extractable amount
    excessive_extraction = extracted_melt >= max_extractable_melt_fraction
    melt_for_crust_generation[excessive_extraction] = 0

    # (2d) Update the remaining melt in the rock
    # remaining_melt = available_melt - melt_for_crust_generation
    meltField.data[:,0] -= melt_for_crust_generation

    # Keep track of how much melt has been extracted in total from a given rock
    cumulativeMelt.data[:,0] += melt_for_crust_generation
    # total_extracted_melt = extracted_melt + melt_for_crust_generation
    
    # overshoot_extraction = total_extracted_melt > max_extractable_melt_fraction
    # total_extracted_melt[overshoot_extraction] = max_extractable_melt_fraction
    overshoot_extraction = cumulativeMelt.data[:,0] > max_extractable_melt_fraction
    cumulativeMelt.data[overshoot_extraction] = max_extractable_melt_fraction

    return melt_for_crust_generation



# %%
def avg_val_at_cell_centres(swarm, swarmVar):
    # Extract data
    material_data = swarmVar.data[:, 0]
    owning_cell_ids = swarm.owningCell.data[:, 0]
    
    # Sort the data by owning cell IDs
    sorted_indices = np.argsort(owning_cell_ids)
    sorted_material_data = material_data[sorted_indices]
    sorted_owning_cell_ids = owning_cell_ids[sorted_indices]
    
    # Find unique owning cell IDs and their start indices
    unique_ids, unique_idx = np.unique(sorted_owning_cell_ids, return_index=True)
    
    # Use np.add.reduceat to sum the groups of sorted_material_data based on unique_idx
    sums = np.add.reduceat(sorted_material_data, unique_idx)
    
    # Calculate counts for each unique ID group
    counts = np.diff(np.append(unique_idx, len(sorted_material_data)))
    
    # Compute the means for each group
    cell_avg = sums / counts
    
    return cell_avg


# %%
def melt_crustal_thickness(cell_melt, dz, resx, resy):
    melt_frac_grid = np.flipud(cell_melt.reshape(resy, resx))
    
    melt_crust_thickness = np.sum(melt_frac_grid, axis=0)*dz

    return melt_crust_thickness


def emplace_crust(CT_f, s_f, min_melt_in_rock, max_extractable_melt_fraction):
    
    ### move all particles downwards for crust emplacement
    move_particles_down(Model.swarm, CT_f)

    # ### move all particles downwards 
    # with Model.swarm.deform_swarm():
    #     Model.swarm.particleCoordinates.data[:,1] -= CT_f(Model.swarm.particleCoordinates.data[:,0])

    ### update material
    Model.materialField.data[:,0][ (Model.swarm.data[:,1] < s_f(Model.swarm.data[:,0])) & 
                                   (CT_f(Model.swarm.particleCoordinates.data[:,0]) > 0) &
                                   (Model.materialField.data[:,0] == air.index )  ] = newCrust.index


    ### update melting region of mantle
    Model.materialField.data[:,0][((Model.materialField.data[:,0] == mantle.index) |
                                   (Model.materialField.data[:,0] == mantleLithosphere.index)) 
                                   & (Model.meltField.data[:,0] > min_melt_in_rock )
                                  ] = partiallyMoltenMantle.index
    
    ### update region that has stopped melting
    Model.materialField.data[:,0][ (Model.materialField.data[:,0] == partiallyMoltenMantle.index) 
                                 & (Model.meltField.data[:,0] <=  min_melt_in_rock ) 
                                 ] = mantle.index
    

    ### update region that has reached the max melt
    Model.materialField.data[:,0][ ((Model.materialField.data[:,0] == mantle.index) |
                                   (Model.materialField.data[:,0] == mantleLithosphere.index) |
                                   (Model.materialField.data[:,0] == partiallyMoltenMantle.index)) 
                                    & (Model.cumulativeMelt.data[:,0] >= max_extractable_melt_fraction ) 
                                 ] = depletedMantle.index
    
def move_particles_down(particles, new_crust_thickness_f):
    ### These particles are only introduced at the end of the extension phase
    try:
        with particles.deform_swarm():
            particles.data[:,1] -= new_crust_thickness_f(particles.data[:,0])

    except:
        pass


    

# %% [markdown]
# ### To make parallel safe:
# - [ ] _cell_meltFrac_ in _generate_new_crust_ needs to be gathered and reshaped into the initial grid
# - [ ] Try to multiply _dz_ in _generate_new_crust_ over the grid in the _z_ direction in case the mesh is deformed

# %%
def generate_new_crust():

    min_melt_in_rock = 0.01
    min_extractable_melt_fraction = 0.05
    max_extractable_melt_fraction = 0.5

    melt_for_crust_generation = update_melt_fractions(Model.meltField, Model.cumulativeMelt,  min_melt_in_rock, min_extractable_melt_fraction, max_extractable_melt_fraction)

    Model.meltExtract.data[:,0] = melt_for_crust_generation

    cell_melt = None

    unique_x = None

    cell_melt_frac = avg_val_at_cell_centres(Model.swarm, Model.meltExtract)
    # cell_melt_frac = Model.meltExtract.evaluate(centroid.data)[:,0]


    cell_material_field = avg_val_at_cell_centres(Model.swarm, Model.materialField)
    # cell_material_field = Model.materialField.evaluate(centroid.data)[:,0]


    ### sort if in parallel
    if uw.mpi.size > 1:
        ''' needs to be sorted on root proc '''
        root_melt_frac = gather_data( cell_melt_frac )
        
        root_coords_x = gather_data( centroid.data[:,0] )
        root_coords_y = gather_data( centroid.data[:,1] )
    
        root_mat = gather_data( cell_material_field )
    
        if uw.mpi.rank == 0:
            ### sort unique x values
            unique_x = np.sort(np.unique( root_coords_x ))
            
            ### gather the melt frac data and coordinates 
            gathered_data = np.column_stack([root_coords_x, root_coords_y, root_melt_frac, root_mat])
            '''sort the data back into the original grid shape'''
            sorted_indices = np.lexsort((gathered_data[:, 0], gathered_data[:, 1]))
            sorted_data = gathered_data[:][sorted_indices]
            melt_extract_sorted = sorted_data[:,2]
    
            mat_sorted = sorted_data[:,3]
    
            cell_melt = melt_extract_sorted
    
            # melt_grid = np.flipud(cell_melt.reshape(Model.elementRes[1], Model.elementRes[0]))
            # mat_grid  = np.flipud(mat_sorted.reshape(Model.elementRes[1], Model.elementRes[0]))
    
            # if Model.step > 0:
            #     import matplotlib.pyplot as plt
            #     fig, ax = plt.subplots(nrows=2, ncols=1)
        
            #     ax[0].imshow(mat_grid, extent=(Model.minCoord[0].m, Model.maxCoord[0].m, Model.minCoord[1].m, Model.maxCoord[1].m))
            #     ax[1].imshow(melt_grid, extent=(Model.minCoord[0].m, Model.maxCoord[0].m, Model.minCoord[1].m, Model.maxCoord[1].m) )
        
            #     plt.savefig(f'./{Model.outputDir}/mat+melt_grid_{Model.checkpointID}.pdf', bbox_inches='tight' )

        
    
        uw.mpi.comm.Barrier()
    
        cell_melt = uw.mpi.comm.bcast(cell_melt, root=0)
        unique_x  = uw.mpi.comm.bcast(unique_x, root=0)

    else:
        cell_melt = cell_melt_frac
        unique_x  = np.unique(centroid.data[:,0])
        


    
    ### need to change in case dz varies
    dz = np.diff(np.unique(centroid.data[:,1]))[0]
    
    
    melt_CT = melt_crustal_thickness(cell_melt, dz, Model.mesh.elementRes[0], Model.mesh.elementRes[1])
    
    x_vals = unique_x


    ## create interp to represent the new crustal thickness
    CT_f = scipy.interpolate.interp1d(x_vals, melt_CT, fill_value="extrapolate")

    ### create interp of surface 
    s_f = scipy.interpolate.interp1d(Model.Surface_tracers.data[:,0], Model.Surface_tracers.data[:,1], fill_value="extrapolate")

    ### emplace the 
    emplace_crust(CT_f, s_f, min_melt_in_rock, max_extractable_melt_fraction)

    try:
        ### Have to manually move tracers down due to melt emplacement at surface 
        move_particles_down(Model.new_crust_tracers, CT_f)
        move_particles_down(Model.sediment_tracers, CT_f)

    except:
        pass

# %%
generate_new_crust()

# %%
if GEO.nProcs == 1 & render:
    Fig = vis.Figure(figsize=(1200, 400))
    Fig.Points(Model.swarm, Model.materialField, fn_size=1.5)
    Fig.show()

# %%
# if GEO.nProcs == 1:
#     lz = m.lithologies.katz.lherzolite()
    
#     p = GEO.dim(pressure, u.pascal).m/1e9
    
#     tliq = np.zeros(np.shape(p))
#     tsol = np.zeros(np.shape(p))
    
#     for i in range(len(p)):
#         tliq[i] = lz.TLiquidus(p[i])
#         tsol[i] = lz.TSolidus(p[i])
    
#     fig = plt.figure(figsize=(8,4))
#     plt.plot(tliq,GEO.dimensionalise(depths, u.kilometer), label="T_s")
#     plt.plot(tsol,GEO.dimensionalise(depths, u.kilometer), label="T_l")
#     plt.plot(GEO.dimensionalise(temperature, u.degC),GEO.dimensionalise(depths, u.kilometer), label="Geotherm")
#     plt.xlabel("Temperature (C)")
#     plt.ylabel("Depth (km)")
#     plt.legend()
#     plt.ylim(110, 0)
#     plt.show()

# %% [markdown]
# ## Solver options

# %%
Model.solver.set_inner_method("superludist")
Model.solver.set_penalty(1e6)
Model.solver.options.scr.ksp_type="cg"


GEO.rcParams["initial.nonlinear.tolerance"] = 1e-2
GEO.rcParams["advection.diffusion.method"] = "SLCN"

# %%
Model.post_solve_functions['advect_PTs']        = advect_tracers_vertically
Model.post_solve_functions['crust_production']  = generate_new_crust
Model.post_solve_functions['surface_processes'] = diffusive_surface

# %% [markdown]
# # Run the Model

# %%
Model.run_for(duration=extension_duration-(2.0001*u.megayear), checkpoint_interval=0.5*u.megayear)



# %%

# %%
### add in tracers over the melt area

# projMaterial = Model.materialField.evaluate(Model.mesh)

# new_crust_coords_x = gather_data(Model.mesh.data[:,0][(projMaterial[:,0] == newCrust.index)], bcast=True)
# new_crust_coords_y = gather_data(Model.mesh.data[:,1][(projMaterial[:,0] == newCrust.index)], bcast=True)

# new_crust_coords = np.column_stack([new_crust_coords_x, new_crust_coords_y])



# %%
#### Extract coordinates of new crust

extractedCoords_x = gather_data(Model.swarm.data[:,0][Model.materialField.data[:,0] == newCrust.index], bcast=True)
extractedCoords_y = gather_data(Model.swarm.data[:,1][Model.materialField.data[:,0] == newCrust.index], bcast=True)

extractedCoords =  np.column_stack([extractedCoords_x, extractedCoords_y])

n_sample = 1500

new_crust_coords = extractedCoords[np.random.randint(0,extractedCoords.shape[0],n_sample)]

Model.add_passive_tracers(name="new_crust", vertices=new_crust_coords, advect=True)


Model.new_crust_tracers.add_tracked_field(Model.pressureField,
                                       name="tracers_press",
                                       units=u.kilobar,
                                       dataType="float")

Model.new_crust_tracers.add_tracked_field(Model.temperature,
                                       name="tracers_temp",
                                       units=u.celsius,
                                       dataType="float")



Model.new_crust_tracers.add_tracked_field(Model.temperature,
                                       name="tracers_time",
                                       units=u.megayear,
                                       dataType="float")

# %%

# %%
#### Extract coordinates of sediment

extractedCoords_x = gather_data(Model.swarm.data[:,0][Model.materialField.data[:,0] == sediment.index], bcast=True)
extractedCoords_y = gather_data(Model.swarm.data[:,1][Model.materialField.data[:,0] == sediment.index], bcast=True)

extractedCoords =  np.column_stack([extractedCoords_x, extractedCoords_y])

n_sample = 500

sediment_coords = extractedCoords[np.random.randint(0,extractedCoords.shape[0],n_sample)]

Model.add_passive_tracers(name="sediment", vertices=sediment_coords, advect=True)



Model.sediment_tracers.add_tracked_field(Model.pressureField,
                                       name="tracers_press",
                                       units=u.kilobar,
                                       dataType="float")

Model.sediment_tracers.add_tracked_field(Model.temperature,
                                       name="tracers_temp",
                                       units=u.celsius,
                                       dataType="float")



Model.sediment_tracers.add_tracked_field(Model.temperature,
                                       name="tracers_time",
                                       units=u.megayear,
                                       dataType="float")


# %%
def save_tracer_PT_data(tracers, tracername):
    pressure_data    = gather_data( Model.pressureField.evaluate(tracers.data) )
    temperature_data = gather_data( Model.temperature.evaluate(tracers.data) )
    time_data        = gather_data( np.repeat(Model.time.m, tracers.data.shape[0]) )  # gather_data( Model.timeField.evaluate(tracers.data) )
    tracer_id_data   = gather_data( tracers.global_index.data )
    if uw.mpi.rank == 0:
        data_to_save = np.column_stack([tracer_id_data, GEO.dim(time_data, u.megayear).m, GEO.dim(temperature_data, u.kelvin).m-273.15, GEO.dim(pressure_data, u.kilobar).m ])


        filePath = f'{Model.outputDir}/PT_data/'
        
        os.makedirs(filePath, exist_ok=True)

        ### Create columns of file if it doesn't exist
        file_name = f'{filePath}PTt-{tracername}-{Model.step}.txt'

        np.savetxt(file_name, data_to_save)

def save_PT_data():
    save_tracer_PT_data(Model.sediment_tracers, 'sediments')
    save_tracer_PT_data(Model.new_crust_tracers, 'newCrust')



Model.post_solve_functions['save_PT_data'] = save_PT_data

# %%
Model.run_for(duration=2.0001*u.megayear, checkpoint_interval=0.5*u.megayear)

# %%
#### Change BC to free slip everywhere (open base still) for 20 Myr

# %%
Model.set_velocityBCs(left=[0., None],
                       right=[0., None],
                       top=[None, 0.], )



# %%
Model.run_for(duration=static_duration + 0.0001*u.megayear, checkpoint_interval=0.5*u.megayear)

# %%
#### Compress the model for 10 Myr

# %%
compression_condition = [(Model.y < GEO.nd(0. * u.kilometre), GEO.nd(-extension_rate)),
                       (True, GEO.nd(-extension_rate) + Model.y * (GEO.nd((2. * extension_rate) / GEO.nd(Model.maxCoord[1]))))]


right_wall_vbc = fn.branching.conditional(compression_condition)


# %%
Model.set_velocityBCs(left=[0., None],
                       right=[right_wall_vbc, None],
                       top=[None, 0.], )


# %%
Model.run_for(duration=collision_duration+0.0001*u.megayear, checkpoint_interval=0.5*u.megayear)

