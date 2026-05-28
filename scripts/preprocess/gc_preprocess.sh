#!/bin/bash
# 
#Script for selecting and moving raw GC data from OutputDir to submission data repository.
#script should be rerun with each inpath/outpath combination to process raw GC data.

module load cdo/1.9.10_oel8
module load nco/5.0.1

#select pair for respective sensitivity run
year=2019
months=$(seq -f "%02g" 1 12)

#inpath=/net/fs03/d0/emroy/GCrundir/run0051/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/NEST

#inpath=/net/fs03/d0/emroy/GCrundir/run0052/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/ANT3-ENA

#inpath=/net/fs03/d0/emroy/GCrundir/run0053/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/ANT10-ENA

#inpath=/net/fs03/d0/emroy/GCrundir/run0054/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/ANT7-ENA

#inpath=/net/fs03/d0/emroy/GCrundir/run0052_NEscale/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/ANT3

#inpath=/net/fs03/d0/emroy/GCrundir/run0053_NEscale/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/ANT10

#inpath=/net/fs03/d0/emroy/GCrundir/run0054_NEscale/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/ANT7

#inpath=/net/fs03/d0/emroy/GCrundir/run0051_FTscale/OutputDir
#outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/CFTSCALMAX

inpath=/net/fs03/d0/emroy/GCrundir/run0051_FTscale_tuned/OutputDir
outpath=/home/emroy/fs03/paper_staging/RegionalBoxModel/data/GCdata/GC_output/CFTSCAL

#First subsetting SpeciesConcHF - need to output SpeciesConcVV_Hg0

cd $outpath

for month in $months
do
  cdo selvar,SpeciesConcVV_Hg0 ${inpath}/GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4 GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4
  cp ${inpath}/GEOSChem.DryDep.${year}${month}01_0000z.nc4 .
  cp ${inpath}/GEOSChem.MercuryEmis.${year}${month}01_0000z.nc4 .
  cp ${inpath}/HEMCO_diagnostics.${year}${month}010000.nc .
  cp ${inpath}/GEOSChem.ProdLoss.${year}${month}01_0000z.nc4 .
done

#create aggregated full year concentration file (useful for supplement)
cdo mergetime GEOSChem.SpeciesConcHF.${year}*0000z.nc4 GEOSChem.SpeciesConcHF.${year}.nc4

