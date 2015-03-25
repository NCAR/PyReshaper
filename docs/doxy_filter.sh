#!/bin/sh
DOXYPYPY=$(which doxypypy)
if [ -x ${DOXYPYPY} ]; then
	$DOXYPYPY -a -c $1
else
	cat $1
fi
