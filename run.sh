#!/bin/bash

ligand_library=$1
box_center=$( echo $2 | sed 's/ //g' )
box_size=$( echo $3 | sed 's/ //g' )
top_n_scores=$4
force_field=$5
docking_method=$6
flexible_sidechains='Empty'

receptor=$(find * -type f | grep -v 'tapisjob')                                                                                   

# set flexible sidechains, if any
#FLEX=${flexible_sidechains}
#if [ -z $FLEX ]; then
#    FLEX='Empty'
#fi

# mpi settings - set VAL based on number of nodes
if [[ "${_tapisNodes}" > 1 ]]; then
    VAL=1
else 
    VAL=0
fi


# Log commands, timing, run job
echo -n "starting: "
date

echo "================================================================"
echo "python3 /autodock-src/autodock.py \\"
echo "    --receptor ${receptor} \\"
echo "    --ligand_library ${ligand_library} \\"
echo "    --center=\"${box_center}\" \\"
echo "    --size \"${box_size}}\" \\"
echo "    --number ${top_n_scores} \\"
echo "    --forcefield ${force_field} \\"
echo "    --docking ${docking_method} \\"
echo "    --sidechains ${flexible_sidechains} "
echo "================================================================"

MV2_ENABLE_AFFINITY=0 MV2_SMP_USE_CMA=${VAL} \
python3 /autodock-src/autodock.py \
     --receptor ${receptor} \
     --ligand_library ${ligand_library} \
     --center="${box_center}" \
     --size "${box_size}" \
     --number ${top_n_scores} \
     --forcefield ${force_field} \
     --docking ${docking_method} \
     --sidechains ${flexible_sidechains}

echo -n "ending: "
date
