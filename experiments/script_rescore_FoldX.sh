#Script to run Foldx for a pair of protein+ligand pdb files
#v2.0
ligID=$1
protID=$2
complex=$3
cname=$4

#cname=${complex%.*}


#Call FoldX
echo -e "\nFoldX will compute the energy between chains $ligID (ligand) and $protID (receptor).\n" 
cp $complex ./${cname}.pdb
./foldx --command=AnalyseComplex --pdb=${cname}.pdb --analyseComplexChains=$ligID,$protID

summary=./Summary_${cname}*.fxout

be=$(sed -n "10p" $summary | awk '{print $6}')
echo -e "\nBinding energy for $lname = $be (Kcal/mol)\n"
echo "$cname, $be" >> FoldX_results.csv
rm ./*.fxout
rm  ./${cname}.pdb
#find . -type f \( -name "*.fxout" -o -name "*_complex.pdb" -o -name "*.pdbqt" \) -delete