# Extrapolation-Boundary-Condition
Extrapolation Boundary Condition is implemented in OpenFOAM. This boundary condition is strongly recommended to use at the outlet of the flow domain. Any flow properties can be extrapolated from the internal field (p U T alpha .. etc). \
Extrapolation boundary condition doesnot take any inputs from the user and hence can be simply implemented as, \
patch_name \
{ \
      type      extrapolation; \
} \
For more details please have a look into into the initial comments on extrapolation/extrapolationFvPatchField.H file.

# Compilation and linking as library
1. create a Make folder outside the extrapolation directory containing "files" and "options" files.
2. In the files, add the path of extrapolationFvPatchFields.C in the first line.
3. Second line goes as, LIB = $(FOAM_USER_LIBBIN)/libextrapolation
4. Run wmake and the last line of the successful compilation looks like, -o /FOAM_USER_DIR/platforms/linux64GccDPInt32Opt/lib/libextrapolation.so

# Linking library through Control Dict
To link this library with the case files of OpenFOAM, Open ControlDict file and add the following line at the end of the file, \
libs ("libextrapolation.so");
