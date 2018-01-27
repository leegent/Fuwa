# Fuwa

A decision-tree-based fast variant caller.

## Install
You can directly download the complete package using the following command:
```bash
wget -c http://rubywu.cn/public/fuwa.gz
```
and then 
```bash
gunzip fuwa.gz
cd fuwa
sh install.sh
```

Alternatively, you can clone this repository to your computer and `make` the executable file.

The reference files and the gziped dbSNP file can be downloaded using the following commands:
```bash
# reference
wget http://rubywu.cn/public/reference.tar.bz2
# dbSNP
wget http://rubywu.cn/public/dbsnp141.gz
```

**Attention:** A 64-bit Linux OS is required for `make`. Type `getconf LONG_BIT` to confirm your OS bits.

## Usage: 
    fuwa [options]
	-i input     input bam file
	-d dbSNP     dbSNP gz file
	-r ref_dir   reference directory
	-o output    [optional] output file name without extension (default: input file name)
	-s SNP_qual  [optional] Filtering threshold of quality score for SNPs (range: [0, 1], default: 0.8)
	-q qual      [optional] variant qual filtering threshold (range: [0, 1], default: 0.6)
	-m           [optional] the sample is male. By default the sample is considered female

## Example
```
./fuwa -i input.bam -o output -r ./reference/ -d ./dbSNP141.gz -q 0.7 -s 0.7
```
