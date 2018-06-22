# Since we run coverage inside docker, we have to fix the paths from the coverage
# file so that they point to the files on the machine instead of inside docker.
cat $1 | sed "s/\/home\/user/workers/g" > "$1.fixed"