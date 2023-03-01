#!/bin/bash
#seting
crdfile="md0.crd"
topfile="pep.top"
startframe=1
endframe=200
interval=1
top=$(pwd)
###################################################

cpptraj > trajout << EOF
parm ${topfile}
trajin ${crdfile} ${startframe} ${endframe} ${interval}
trajout md.crd crd nobox
EOF

cat > energy_rms.in << EOF
&energy_rms_cntrl

snapshot_start = 1
snapshot_end = 200
snapshot_step =1         !number of snapshot each step
snapshot_out = 20         !number of output snapshot

natom_start_pro = 1      !default value 1
natom_end_pro   = 14476
natom_start_lig = 14477
natom_end_lig = 14531

totenergyfile = 'totenergy.dat'
energyfile = 'energy.dat'

ifbox = 0        !the box in crdfile, 0:NO; 1:YES default value 1
cut = 999        !default value 999
dielc= 1.0       !default value 1
&end
EOF

gfortran energy_rms.f90 -o energy_rms
./energy_rms -i energy_rms.in -o energy_rms.out -p ${topfile} -c md.crd

rm energy_rms energy_rms.in md.crd trajout


