echo "#####  start flow analysis..."

for pid in {0..2}
do
echo "#####  starting pid " $pid
#root -l "runNewJuly2016.C( $pid )"
root -l "runNewOct2017.C( $pid )"
done
