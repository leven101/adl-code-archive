ss=`wc -l wcoutput`
for word in $ss
do 
  lines=$word
  break;
done
echo "Total lines in file: $lines"
src=-1;
testPts=36;
for (( i=$testPts; i<=$lines; i+=$testPts)) 
do
  if [ $src -eq 33 ]
  then
    src=-1
  fi
  echo "Source $((src+=1))"
  head -n $i wcOutput | tail -n $testPts | awk '{print $6-1"\t"$4}' | ../perf.src/perf -acc -prf -roc
done
