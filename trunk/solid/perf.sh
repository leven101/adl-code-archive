ss=`wc -l wcoutput`
for word in $ss
do 
  lines=$word
  break;
done
ss=`tail -n 1 wcoutput`
i=0;
for word in $ss
do
  i=$((i+1));
  if [ $i -eq 1 ]
  then 
    wkcls=$word
  fi
  if [ $i -eq 2 ]
  then
    totEps=$word
    break
  fi
done
src=-1;
testPts=$1;
echo "Total lines in file: $lines"
echo "Total epochs: $totEps. Testing epochs: $testPts"
for (( i=$totEps; i<=$lines; i+=$totEps)) 
do
  echo "Source $((src+=1))"
  head -n $i wcOutput | tail -n $testPts | awk '{print $6-1"\t"$4}' | ../perf.src/perf -acc -prf -roc
done
