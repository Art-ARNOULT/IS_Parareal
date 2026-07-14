for i in $(seq 1 1 80);
do

    # echo $i
    ./exe_fine $i &

done
