#!/bin/bash
set -x

cd $WEST_STRUCT_DATA_REF
python $WEST_SIM_ROOT/common_files/init_pcoords.py
paste $WEST_STRUCT_DATA_REF/dist1.dat $WEST_STRUCT_DATA_REF/dist2.dat $WEST_STRUCT_DATA_REF/angle1.dat > $WEST_PCOORD_RETURN
