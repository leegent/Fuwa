## Usage: 
    fuwacall [options]
	-i input     input bam file
	-d dbSNP     dbSNP gz file
	-r ref_dir   reference directory
	-o output    [optional] output file name without extension (default: input file name)
	-q qual      [optional] variant qual filtering threshold (range: [0, 1], default: 0.8)
	-m           [optional] the sample is male. By default the sample is considered female

## example:
```
fuwacall -i input.bam -o output -r reference/ -d dbSNP141.gz -m
```

The reference dicretory and gziped dbSNP can be downloaded [here](http://pan.baidu.com/s/1pKBRIPl)