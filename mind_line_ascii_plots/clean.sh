#!/bin/bash
# Cleans file directory tree of compiled files

# Remove output dir(s)
cd `dirname $0`

rm -rf build

# remove data directory
rm -rf data

# Let the User know the directory is clean
echo "The project is clean and ready to build!"
