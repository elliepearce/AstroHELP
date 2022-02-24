# AstroHELP

This interface takes MCFOST output and runs it through ASTROCHEM, a code which computes the abundances of species in an interstellar medium as a function of time,  to provide better models of accretion and protoplanetary discs. The output from this interface can then be run again through MCFOST to further improve models. The main aim of this interface is to improve efficency of computing more complex models that can be used for more accurate prediction of abundances of within interstellar mediums.

### Required Input Format:
The MCFOST output required to run AstroHELP is the following:
- UV Field
- Column density
- Gas density
- Temperature
- Dust mass desnity
- Grain sizes
- Dust particle density
- Grid
#### These files must all be in the same directory.

(To generate UV field, include line: -output_UV_field. To generate column density, include line: -cd)
