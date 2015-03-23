fin="ibccOut"  
ss=`wc -l $fin`
for word in $ss
do 
  lines=$word
  break;
done
ss=`tail -n 1 $fin`
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
    totEps=$((word+1))
    break
  fi
done
src=-1;
testPts=$1;
c2c=$2
echo "Total lines in file: $lines"
echo "Total epochs: $totEps. Testing epochs: $testPts"
for (( c=$totEps; c<=$lines; c+=$totEps ))
do
   echo "Welcome $c times"
done
for (( i=$totEps; i<=$lines; i+=$totEps ))
do
   echo "Source $((src+=1))"
   head -n $i $fin | tail -n $testPts | egrep "(\t$c2c\t[0-9]$)|(\t$c2c$)" \
    | awk '{print(c2c != $NF ? 0 : 1)"\t"$(c2c+3)}' c2c=$2 | ../perf.src/perf -acc -prf -pre -rec -roc 
done
