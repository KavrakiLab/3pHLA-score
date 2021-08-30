#!/bin/bash

#Script to assess the binding energy of a complex using GradDock
#v1

ligand=$1
lname=${ligand%.*}
receptor=$2
rname=${receptor%.*}
directory=$3

echo -e "________________________________________________________________________________"
echo -e "\nEvaluate energy with GradDock"
echo -e "________________________________________________________________________________"
echo -e "Ligand: $1"
echo -e "Receptor: $2"

echo -e "Lname: $lname"
echo -e "Rname: $rname"
echo -e "file looking $directory/$ligand"

#Preprocessing the ligand - adding all hydrogens
echo -e "\nConverting ligand to pdbqt format - adding all hydrogens..."
/data/reduce.3.23.130521 $directory/$ligand > ${lname}_h.pdb

#Identifying alfa chain and delete beta chain
protID=$(grep ATOM $directory/$receptor | sed -n 1p | awk '{print $5}')
awk '($5 == id)' id="$protID" $directory/$receptor > ${rname}_fixed.pdb

#Removing all ions and additional molecules from the pdb
sed -i '/HETATM/d' ${rname}_fixed.pdb

#Preprocessing the receptor - adding all hydrogens
echo -e "\nConverting receptor to pdbqt format - adding all hydrogens..."
/data/reduce.3.23.130521 ${rname}_fixed.pdb > ${rname}_h.pdb

source ./build  &&

#Format the ligand in the way GradDock score requires (first line echo name, last line ENDMDL)
echo echo ${lname}_h.pdb > ${lname}_formatted.pdb
cat ${lname}_h.pdb | sed \$d | sed '/^ATOM/!d' >> ${lname}_formatted.pdb
echo "ENDMDL" >> ${lname}_formatted.pdb
#Format empty file (for scoring receptor alone) in the way GradDock score requires
echo -e "echo empty\nENDMDL" > empty.pdb

echo -e "\nScoring the complex with GradDock..."
./score D47 ${rname}_h.pdb < ${lname}_formatted.pdb > ${rname}_complex_results.txt

res1=$(sed -n -e '/^ Total weighted score:  /p' ${rname}_complex_results.txt | awk '{print $4}')
echo -e "\nComplex score: $res1"

echo -e "\nScoring the receptor with GradDock..."
./score D47 ${rname}_h.pdb < empty.pdb > ${rname}_results.txt

res2=$(sed -n -e '/^ Total weighted score:  /p' ${rname}_results.txt | awk '{print $4}')
echo -e "\nReceptor score: $res2"

echo -e "\nScoring the ligand with GradDock..."
./score D47 ${lname}_h.pdb < empty.pdb > ${lname}_results.txt

res3=$(sed -n -e '/^ Total weighted score:  /p' ${lname}_results.txt | awk '{print $4}')
echo -e "\nLigand score: $res3"

r=$(echo "($res1-$res2-($res3))"| bc -l)
echo -e "\nBinding energy: $r"

echo -e "$directory/$receptor, $res1, $r" >> GD_results.csv

rm ${lname}_h.pdb
rm ${rname}_fixed.pdb
rm ${rname}_h.pdb
rm ${lname}_formatted.pdb
rm ${rname}_complex_results.txt
rm ${rname}_results.txt
rm ${lname}_results.txt

echo -e "________________________________________________________________________________"
