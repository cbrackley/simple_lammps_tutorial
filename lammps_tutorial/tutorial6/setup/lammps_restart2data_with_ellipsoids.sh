#!/bin/bash

path2lammps="~/work/lmp-7Dec15_openmpi_SL7"

# In the current version of LAMMPS, the write_data command does not 
# output the ellipsoids

if [[ ! $# == 2 ]]; then
    echo "Usage : lammps_restart2data_with_ellipsoids.sh input.restart output.data"
    exit 1
fi

infile=$1
outfile=$2

if [[ ! -f $infile ]];then
    echo "ERROR : Cannot find file $infile"
    exit 1
fi

if [[ -f $outfile ]]; then
    echo "ERROR : File $outfile already exists"
    exit 1
fi

#####
#####

## Notes : 
# 1. LAMMPS complains that angle coeff must be set. 
#    The script below assumes there are two angle 
#    types - need to edit this if there are more.
#    Since no simulation is run, and coeff and not 
#    output to the file, doesn't matter what these are.
# 2. The script also assumes that all atoms are ellipsoids.
#    Not sure what would happen if we attempt to compute 
#    quaternions for atoms which do not have that property.
# 3. The compute shape seems to give the radius rather than 
#    the diameter as stated in the LAMMPS manual. The data
#    file must specify radius.


eval $path2lammps << EOF

atom_style hybrid angle ellipsoid
read_restart $infile

angle_style   cosine
angle_coeff   1 20.0
angle_coeff   2 20.0

compute shap all property/atom shapex shapey shapez
compute quat all property/atom quatw quati quatj quatk
dump d2 all custom 1 temp_ell id c_shap[1] c_shap[2] c_shap[3] c_quat[1] c_quat[2] c_quat[3] c_quat[4]
dump_modify d2 first yes
run 0

write_data temp_data nocoeff

EOF


### Now count the number of ellipsoids in the data file
Nell=$( awk 'BEGIN{A=0;V=0;}{
if ($1=="Velocities") {A=0;V=1;}
if (A==1&&$1!="") {c+=$7}
if ($1=="Atoms") {A=1;}
}END{print c}' temp_data )

### And make the output file
awk -v c=$Nell '{
print;
if ($2=="atoms") {print c " ellipsoids";}
}' temp_data > $outfile

echo "" >> $outfile
echo " Ellipsoids" >> $outfile
echo "" >> $outfile

awk '{if (NR>9) print $1,$2*2,$3*2,$4*2,$5,$6,$7,$8}' temp_ell >> $outfile
# radii are multipled by 2

## clean up
rm temp_data temp_ell
