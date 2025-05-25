if [ "$#" -ne 2 ]; then
	echo "You need to specify the test and amount of processes and in that order."
	exit 1
fi

echo "mpic++ $1 -o _test"
echo "mpiexec -n $2 _test"
echo "rm _test"

mpic++ $1 -o _test
mpiexec -n $2 _test
rm _test