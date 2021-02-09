#!/bin/bash
#=======================================================================
#  Crude parallelization of RFIClean for processing filterbank files.
#
#           Yogesh Maan, Dec. 2019
#=======================================================================
#


#------------------------------------------------------------------------------
if [ $# -lt 4 ] ; then
 echo "========================================================================================"
 echo " Usage: "
 echo " crp_rficlean_fil.sh  <out_filename>  <flag-file>  <N_parallel>  <in_filename>" 
 echo " "
 echo "========================================================================================"
 exit 0
fi
outfile=$1
flagfile=$2
Np=$3
infile=$4
#------------------------------------------------------------------------------
rtag=${outfile%*.fil}  # tag for part-files

##-------------------------------
## read flags from the flag-file
block=`sed -n 1p $flagfile`
flag=`sed -n 2p $flagfile`
#echo "block-size: $block"
#echo "flag:   '$flag'"
##-------------------------------


# deduce number of blocks to be processed in each part
nsamples=`header $infile | grep "Number of samples" | awk '{print $5}'`
Nb=`echo "scale=0; $nsamples/($Np*$block) + 1" | bc`
#echo "nsamples, Nblocks: $nsamples, $Nb, $Np, $block ---"



echo "Starting parallel-RFIClean processing at: "
date
echo ""

##-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# prepare and execute the commands for all parts
mkdir -p RFIClean_ps
combine_cmd="cat "
remove_cmd="rm  "
count=0
for (( c=1; c<=$Np; c++ )); do
  bst=$(($count*$Nb + 1))
  count=$((count+1))
  #echo $count
  tempout=${rtag}_part${count}.fil
  psfile=${rtag}_part${count}.ps
  #-------------------------------------
  if [ $count -eq 1 ]; then
    cmd="rficlean -t $block $flag -o  $outfile -ps RFIClean_ps/$psfile  $infile -bst $bst -nbl $Nb &"
  elif [ $count -eq $Np ]; then
    combine_cmd="${combine_cmd} ${tempout}"
    remove_cmd="${remove_cmd} ${tempout}"
    cmd="rficlean -t $block $flag -o  $tempout -ps RFIClean_ps/$psfile  $infile -bst $bst -headerless &"
  else
    combine_cmd="${combine_cmd} ${tempout}"
    remove_cmd="${remove_cmd} ${tempout}"
    cmd="rficlean -t $block $flag -o  $tempout -ps RFIClean_ps/$psfile  $infile -bst $bst -nbl $Nb -headerless &"
  fi
  echo $cmd
  eval $cmd
done
wait
##-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


# Now combine all parts and then remove the parts
combine_cmd="${combine_cmd} >> ${outfile}"
#echo $combine_cmd
#echo $remove_cmd
cmd="$combine_cmd && $remove_cmd"
echo $cmd
eval $cmd


echo "Finished parallel-RFIClean processing at: "
date
echo ""

exit 0

