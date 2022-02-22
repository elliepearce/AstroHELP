# AstroHELP

This interface takes MCFOST output and runs it through ASTROCHEM inorder to provide better models of accretion and protoplanetary discs. The output from this interface can then be run again through MCFOST to further improve models.

### Required Input Format:
The MCFOST output required is:
- UV Field
- Column density
- Gas density
- Temperature
- Dust mass desnity
- Grain sizes
- Dust particle density
- Grid
#### These files must all be in the same folder/directory.

To generate UV field, include line: -output_UV_field.
To generate column density, include line: -cd
