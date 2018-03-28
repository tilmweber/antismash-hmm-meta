#!/bin/sh
for fn in *.hmm; do
	dn=${fn/.hmm/}
	echo Creating $dn
	mkdir $dn
	mv $fn $dn
done
	