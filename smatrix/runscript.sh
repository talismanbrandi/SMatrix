#!/bin/bash

mkdir -p RUNS
cd RUNS
mkdir -p logs

for n1r in 0 1
do
    for n2r in 0 1
    do
        for n3r in 0 1
        do
            for n4r in 0 1
            do
                for n5r in 0 1
                do
                    for n1i in 0 1
                    do
                        for n2i in 0 1
                        do
                            for n3i in 0 1
                            do
                                for n4i in 0 1
                                do
                                    for n5i in 0 1
                                    do
                                        for nd1 in 0 1
                                        do
                                            n=$((2#${n1r}${n2r}${n3r}${n4r}${n5r}${n1i}${n2i}${n3i}${n4i}${n5i}${nd1}))
                                            mkdir -p configuration_${n}
                                            cd configuration_${n}
                                            cp -r ../../config .
                                            sed -i -e "s/N1R/${n1r}/" config/myModel.conf
                                            sed -i -e "s/N2R/${n2r}/" config/myModel.conf
                                            sed -i -e "s/N3R/${n3r}/" config/myModel.conf
                                            sed -i -e "s/N4R/${n4r}/" config/myModel.conf
                                            sed -i -e "s/N5R/${n5r}/" config/myModel.conf
                                            sed -i -e "s/N1I/${n1i}/" config/myModel.conf
                                            sed -i -e "s/N2I/${n2i}/" config/myModel.conf
                                            sed -i -e "s/N3I/${n3i}/" config/myModel.conf
                                            sed -i -e "s/N4I/${n4i}/" config/myModel.conf
                                            sed -i -e "s/N5I/${n5i}/" config/myModel.conf
                                            sed -i -e "s/ND1/${nd1}/" config/myModel.conf
                                            name=`printf "%.4d" ${n}`
                                            ln -s /home/apaul/beegfs/SMatrix/smatrix/analysis .
                                            cd ../
                                        done
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
cp ../sbatch.sh .
sed -i -e "s/NAME/RUNS-${name}/" sbatch.sh
sbatch sbatch.sh
cd ..