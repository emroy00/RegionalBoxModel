#!/bin/bash
#
#Script for placing all Hg0 concentration data in single file. Useful for regional box model paper.

module load cdo/1.9.10_oel8
module load nco/5.0.1

year=2019
months=$(seq -f "%02g" 1 2)

run=run0045
path=/net/fs03/d0/emroy/GCrundir/${run}/OutputDir/
outpath=${path}/analysis/

cd $path

for month in $months
do
  #ncks -v SpeciesConcVV_Hg0,hyai,hybi,ilev GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4 ${outpath}GEOSChem.SpeciesConcHF.${year}${month}_temp.nc4
  #cdo delvar,AREA GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4 ${outpath}GEOSChem.SpeciesConcHF.${year}${month}_temp.nc4
  cdo selvar,SpeciesConcVV_Hg0 GEOSChem.SpeciesConcHF.${year}${month}01_0000z.nc4 ${outpath}GEOSChem.SpeciesConcHF.${year}${month}_temp.nc4
  #cdo selvar,U,V MERRA2.201706${day}.A3dyn.2x25.nc4 ${outpath}MERRA2.${year}${month}${day}.UV.temp.2x25.nc4
done

cd $outpath
cdo mergetime GEOSChem.SpeciesConcHF.${year}*_temp.nc4 GEOSChem.SpeciesConcHF.${year}.nc4
#cdo mergetime MERRA2.${year}${month}*UV.temp.2x25.nc4 MERRA2.${year}${month}.UV.2x25.nc4
rm GEOSChem.*temp.nc4

