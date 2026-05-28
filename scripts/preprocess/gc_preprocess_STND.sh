#!/bin/bash
# 
#Script for selecting and moving raw GC data from OutputDir to submission data repository.

module load cdo/1.9.10_oel8
module load nco/5.0.1

#select pair for respective sensitivity run
year=2019
months=$(seq -f "%02g" 1 12)

inpath=/net/fs03/d0/emroy/GCrundir/run0045/OutputDir
outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/STND

#First subsetting SpeciesConcHF - need to output SpeciesConcVV_Hg0

cd $outpath

for month in $months
do
  ncks -v SpeciesConcVV_Hg0,hyai,hybi,ilev ${inpath}/GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4 GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4
  #cdo selvar,SpeciesConcVV_Hg0 ${inpath}/GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4 GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4
  #cp ${inpath}/GEOSChem.DryDep.${year}${month}01_0000z.nc4 .
  #cp ${inpath}/GEOSChem.MercuryEmis.${year}${month}01_0000z.nc4 .
  #cp ${inpath}/HEMCO_diagnostics.${year}${month}010000.nc .
done

#create aggregated full year concentration file (useful for supplement)
cdo mergetime GEOSChem.SpeciesConcHF.${year}*0000z.nc4 GEOSChem.SpeciesConcHF.${year}.nc4

