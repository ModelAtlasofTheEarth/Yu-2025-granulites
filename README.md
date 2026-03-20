# New [M@TE](https://mate.science/)! model: 
 _we have provided a summary of your model as a starting point for the README, feel free to edit_
## Section 1: Summary of your model   

**Model Submitter:**  

Ben Steven Knight ([0000-0001-7919-2575](https://orcid.org/0000-0001-7919-2575))

**Model Creator(s):**  

- Ben Steven Knight ([0000-0001-7919-2575](https://orcid.org/0000-0001-7919-2575))  
  
**Model slug:**  

`Yu-2025-granulites` 

(this will be the name of the model repository when created) 

**Model name:**  

_Rift foundering and the generation of non-orogenic granulites during the Mesoproterozoic_  

**License:**  

[Creative Commons Attribution 4.0 International]( https://creativecommons.org/licenses/by/4.0/legalcode.txt)

**Model Category:**  

- model published in study   
  
**Model Status:**  

- completed   
  
**Associated Publication title:**  

_[Rift foundering and the generation of non-orogenic granulites during the Mesoproterozoic](https://doi.org/10.1016/j.epsl.2025.119681)_ 

**Short description:**  

This model was developed to test the generation of non-orogenic granulites due to rifting. We implemented a melt generation and emplacement model and compared with analytical results obtained from the Fraser zone, SW WA.

**Abstract:**  

Mesoproterozoic orogens are unusual in that they commonly preserve a record of high geothermal gradients, low crustal thickness, and limited topography. One such terrane is the Fraser Zone in the Albany–Fraser Orogen (AFO), Western Australia, where granulite-facies rocks record a counterclockwise pressure–temperature (CCW $P–T$) evolution, the drivers of which remain the subject of debate. To understand how the Fraser Zone reached high temperatures followed by burial, a new geochronological and petrological dataset from a metamorphosed gabbronorite was collected. This data places direct temporal constraints on the up-pressure section of the CCW $P–T$ evolution. The gabbronorite preserves an original cumulate texture that crystallised at $\sim 6\ \text{kbar}$ and $\sim 810^{\circ}\text{C}$ that partially recrystallised at conditions of $\sim 9.5\ \text{kbar}$ and $850{-}950^{\circ}\text{C}$, as constrained by conventional thermobarometers and pseudosection modelling. Zircon grains with distinct textural and geochemical characteristics constrain the timing of emplacement, melt crystallisation at $1288 \pm 7\ \text{Ma}$, and subsequent up-pressure recrystallisation at $1284 \pm 7\ \text{Ma}$. These combined $P–T$ and geochronological constraints define a rapid burial path to $\sim 12\ \text{km}$ within $c.\ 4\ \text{Myr}$ of initial crystallisation, requiring a re-evaluation of the previous models involving collision and thickening as the timing of the up-pressure excursion pre-dates the established timing of collision between the Western Australian and South Australian Cratons. Thermomechanical geodynamic modelling elucidates a viable tectonic setting for the generation of the granulites of the Fraser Zone, involving rift foundering triggered by asymmetric extension in the backarc that was terminated by subsequent arc advance. Globally, a similar mechanism may have resulted in the ubiquitous high thermobaric ratio metamorphism, low crustal thickness, and limited elevation, associated with the Mesoproterozoic metamorphic record associated with the assembly of Rodinia.

**Scientific Keywords:**  

- granulites   
- melt   
- emplacement   
  
**Funder(s):**  
- ARC (https://ror.org/05mmh0f86)  
  
## Section 2: your model code, output data  

**No embargo on model contents requested** 

**Include model code:**   

True 

**Model code notes:**   

Model is setup with a python script
Models require underworld2 to be run 

**Include model output data:**   

True 

**Model output data notes:**   

output is primarily h5 files, with 139 timesteps saved. ~8.18 GB of data. 

## Section 3: software framework and compute details   
**Software Framework DOI/URL:**  

Found software: _[underworld2](https://zenodo.org/records/6820562)_ 

**Software Repository:**   

https://github.com/underworldcode/underworld2 

**Name of primary software framework:**  

underworld2 

**Software & algorithm keywords:**  

- python   
- finite element   
- particle in cell   
  
## Section 4: web material (for mate.science)   
**Landing page image:**  

Filename: [graphics/Model_evolution.pdf](https://github.com/user-attachments/files/23468159/Model_evolution.pdf)  
Caption:  Evolution of model at selected timesteps, showing the melt generation and emplacement.  
  
**Animation:**  

Filename: [graphics/animation](https://github.com/user-attachments/assets/ab1c6548-fab2-42ad-bd5c-b0d78f670bb7)  
Caption:  Model evolution showing the generation and emplacement of melt during rifting.  
  
**Graphic abstract:**  

Filename: [None]()  
  
**Model setup figure:**  

Filename: [graphics/Model_setup.pdf](https://github.com/user-attachments/files/23468175/Model_setup.pdf)  
Caption:  Model setup, showing the initial geotherm and material distribution.  
Description:  The 2D model is designed to simulate extension and the emplacement of melt in the crust. The model has a length (x) of 660 km and a height (y) of 140 km. The grid is uniformly spaced at 330 x 70 nodes, producing a grid resolution of 2 km, with 30 particles per cell to track material properties. The model is layered, with a 20 km thick upper crust and 20 km thick lower crust, 80 km thick lithospheric mantle and 10 km thick asthenosphere. A Gaussian plastic strain distribution is initially prescribed across the crust and lithospheric mantle localises deformation and promotes the thinning of the crust during extension. A constant temperature (T = 20 °C) is applied to the top boundary, with no heat flux across the side walls. A Moho temperature of 700 °C is prescribed at a depth of 40 km, which results in a geotherm 17 °C/km across the crust. The lithosphere-asthenosphere temperature is 1375 °C at a depth of 120 km, which is a geothermal gradient of 8.4375 °C/km across the lithosphere. In the asthenosphere, a 0.4 °C/km adiabatic gradient is prescribed, resulting in a temperature of 1380 °C at the bottom of the domain.

  
