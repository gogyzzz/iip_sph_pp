#!/bin/sh

for arg in "$@"
do
  echo "$arg"

sed -i 's/cxmul/CXMUL/g' $arg
sed -i 's/cxadd/CXADD/g' $arg
sed -i 's/cxadd_mul/CXADD_MUL/g' $arg
sed -i 's/cxemul/CXEMUL/g' $arg
sed -i 's/cxediv/CXEDIV/g' $arg
sed -i 's/cxeadd/CXEADD/g' $arg

done
