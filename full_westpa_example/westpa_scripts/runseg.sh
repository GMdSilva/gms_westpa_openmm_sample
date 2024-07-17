#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
  set -x
  env | sort
fi

cd $WEST_SIM_ROOT
mkdir -pv $WEST_CURRENT_SEG_DATA_REF
cd $WEST_CURRENT_SEG_DATA_REF

ln -sv $WEST_SIM_ROOT/common_files/abl1_drude_config.json .

if [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_CONTINUES" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/prod_drude_westpa.py > prod.py
  cp $WEST_PARENT_DATA_REF/seg.rst ./parent.rst
  cp $WEST_PARENT_DATA_REF/seg.pdb ./parent.pdb
elif [ "$WEST_CURRENT_SEG_INITPOINT_TYPE" = "SEG_INITPOINT_NEWTRAJ" ]; then
  sed "s/RAND/$WEST_RAND16/g" $WEST_SIM_ROOT/common_files/prod_drude_westpa.py > prod.py
  cp $WEST_SIM_ROOT/common_files/parent.rst .
  cp $WEST_PARENT_DATA_REF/basis.pdb ./parent.pdb
fi

# Run the dynamics with OpenMM
python prod.py

#Calculate pcoord with MDTraj
python $WEST_SIM_ROOT/common_files/get_pcoords.py

cat dih1.dat > $WEST_DIH1_RETURN
cat dih2.dat > $WEST_DIH2_RETURN

paste dist.dat dist2.dat angle1.dat > $WEST_PCOORD_RETURN
# Clean up
rm -f prod.py seg.dcd dih1.dat dih2.dat dist.dat dist2.dat angle1.dat
