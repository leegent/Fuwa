#!/bin/bash
echo "Unpackage samtools.."
unzip samtools-0.1.19.zip
echo "Make executable file.."
make
chmod a+x ./fuwa
echo "Unpackage reference directory.."
tar -jxvf reference.tar.bz2
echo "Done."