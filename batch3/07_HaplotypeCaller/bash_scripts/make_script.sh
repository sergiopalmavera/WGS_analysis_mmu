# Take the example script and change it according to the desired batch of samples


FROM=$1
TO=$2

script_nm=HaplotypeCaller_xx_to_yy.sh

fl_nm=$(echo $script_nm | sed "s/xx/$FROM/;s/yy/$TO/")

sed -z "s/FROM/FROM=$FROM/ ; s/TO/TO=$TO/" $script_nm > $fl_nm # -z to match only first instance

chmod +x ./$fl_nm

nohup ./$fl_nm &> ./${fl_nm/.sh/.out} &



