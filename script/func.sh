#!/bin/sh

awk '{
split($2,a,"(");
print a[1];
}' $1

