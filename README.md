# Spray Dry with spray_py !
Physics-Based Spray Dryer model, to help tune parameters and guide spray dryer design
Started from a droplet heat model I made for Chad Curtis' University of Washington Chemical Engineering Graduate Mathematics class.


## Model Inputs:
There will be several functionalities of the model, so that it can take either experimental empirical data to generate a matching performance model, or else to generate a predictive model to guide design from first principles. 

Similarly, there are several use-cases and purposes that will process and use data for various inputs/outputs, such as guiding overall spray dryer design and purchasing, or improving slurry formulation to acheive a desired outcome, or even tuning control parameters to design a process control and automation setup. Here are just a list of possibly important design inputs needed to make a physics model of the system:

### Spray Dryer Design Inputs:
 * Dimensions (Width, Cylinder Height, Cone Height)
 * Flow profile (Cocurrent, Countercurrent, Fountain...)
 * Spray Input Location
 * Spray Type and Cone-Shape / velocity calculation profile
 * Collection Mechanism (Cyclone Design?) and location...
 * Gas flow design and baffles/controls 
 * Temperature control mechanism and feedback loop timing (delay parameters)
 * Thermal insulation (across various locations)
 
### Slurry + Material Inputs:
 * Solids Content %wt
 * Solid materials, or else the bulk parameters of the finished material
 * Solvent type (Water, Ethanol...)
 * Slurry viscosity + thinning properties
 * Material thermal and chemical stability (sugars etc to react?)
 * Final PSD Desired (?)
 * Sub-particle size and type (?)
 * dry material "stickiness?"
 * pH (?) 
 
### Process Inputs:
 * Gas content (Air, N2...)
 * Gas flow rate and control mechanism
 * Internal gas pressure (Pos Pressure, Neg Pressure flow)
 * Inlet heating source/type
 * Spraying Pump Type + Flow Rate control
 * Spray Nozzle Setup + Process parameters (RPM for spinning disc)
 
