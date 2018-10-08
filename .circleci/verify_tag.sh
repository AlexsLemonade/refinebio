#!/bin/bash

# This script verifies that the tag triggering this deploy was signed
# by a trusted member of the CCDL.

dongbo_key="138F2211AC85BFB74DB4F59405153C3530E360A7"
rich_key="65EC6E219A655CA8CAC3096890D90E5F0A053C53"
casey_key="DFAA02F5552C553B00CC3DCC31D330047976BAA1"
kurt_key="B539C2CA0D4424660876A9381E1C8D1C2A663250"

trusted_keys="$dongbo_key $rich_key $casey_key $kurt_key"

for key in $trusted_keys; do
    gpg --keyserver pgp.mit.edu --recv-keys $key
done

# If it is not a good key then the exit code is 1, which will cause
# the deploy to fail.
git tag --verify $CIRCLE_TAG
