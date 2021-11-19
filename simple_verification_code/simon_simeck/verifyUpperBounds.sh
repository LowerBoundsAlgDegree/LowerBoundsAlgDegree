#!/bin/zsh

TIMELIMIT=500

echo "Start upper bound computation"
./upperBounds 16 simon $TIMELIMIT > upperbounds_simon32.txt &
./upperBounds 24 simon $TIMELIMIT > upperbounds_simon48.txt &
./upperBounds 32 simon $TIMELIMIT > upperbounds_simon64.txt &
./upperBounds 48 simon $TIMELIMIT > upperbounds_simon96.txt &
./upperBounds 64 simon $TIMELIMIT > upperbounds_simon128.txt &

./upperBounds 16 simeck $TIMELIMIT > upperbounds_simeck32.txt &
./upperBounds 24 simeck $TIMELIMIT > upperbounds_simeck48.txt &
./upperBounds 32 simeck $TIMELIMIT > upperbounds_simeck64.txt &

wait
echo "upper bound computation finished"