#!/bin/bash
#=======================================================================
#  Crude parallelization of RFIClean for processing GMRT data files.
#
#           Yogesh Maan, Dec. 2019
#=======================================================================
#



#------------------------------------------------------------------------------
if [ $# -lt 6 ] ; then
 echo "========================================================================================"
 echo " Usage: "
 echo " crp_rficlean_gm.sh  <out_filename>  <flag-file>  <N_parallel>  <in_filename>  <gm_info>  <extra-flag>" 
 echo " "
 echo " Note: 'rficlean' should be in your PATH."
 echo "========================================================================================"
 exit 0
fi
outfile=$1
flagfile=$2
Np=$3
infile=$4
gminfo_file=$5
extra_flag=$6
#------------------------------------------------------------------------------

rtag=${outfile%*.fil}  # tag for part-files

##-------------------------------
## read flags from the flag-file
block=`sed -n 1p $flagfile`
flag=`sed -n 2p $flagfile`
#echo "block-size: $block"
#echo "flag:   '$flag'"
## read nchan from the gm-info-file
nchans=`sed -n 4p $gminfo_file`
##-------------------------------


# deduce number of blocks to be processed in each part
fsize=$(stat -cL%s "$infile")
nsamples=`echo "scale=0; $fsize/(2*$nchans)" | bc`    ## gmrt data --> 2 bytes/sample
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
  pdffile=${rtag}_part${count}.pdf
  #-------------------------------------
  if [ $count -eq 1 ]; then
    cmd="rficlean -t $block $flag -o  $outfile -ps RFIClean_ps/$psfile  $infile -bst $bst -nbl $Nb -gm $gminfo_file $extra_flag  && ps2pdf  RFIClean_ps/$psfile RFIClean_ps/$pdffile  &&  rm RFIClean_ps/$psfile &"
  elif [ $count -eq $Np ]; then
    combine_cmd="${combine_cmd} ${tempout}"
    remove_cmd="${remove_cmd} ${tempout}"
    cmd="rficlean -t $block $flag -o  $tempout -ps RFIClean_ps/$psfile  $infile -bst $bst -headerless -gm $gminfo_file $extra_flag  && ps2pdf  RFIClean_ps/$psfile RFIClean_ps/$pdffile  &&  rm RFIClean_ps/$psfile &"
  else
    combine_cmd="${combine_cmd} ${tempout}"
    remove_cmd="${remove_cmd} ${tempout}"
    cmd="rficlean -t $block $flag -o  $tempout -ps RFIClean_ps/$psfile  $infile -bst $bst -nbl $Nb -headerless -gm $gminfo_file $extra_flag  && ps2pdf  RFIClean_ps/$psfile RFIClean_ps/$pdffile  &&  rm RFIClean_ps/$psfile &"
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

