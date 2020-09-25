# Ascat-cnv

## Introduction
Prepare ASCAT input files from alleleCount on WES data and run ascat
## Useage
```
Description     
Options
    -dir <s>: dir of bam files
	-outdir <s>: pathway of outdir
	-info <s>: pair information of tumor-normal, for example "normal tumor patiant",first cloumn is normal,second column is tumor,third column is patient
	-alleleCounter <s>: path of alleleCounter
	-multi: run ASCAT with multi-sample segmentation (when shared breakpoints are expected)
	-cpu <i>: cpu used
	-node <s>: nodes used
        -h|?|help : Show this help
```
### HaHa,this pipiline just for me.
### Ofcourse,Other users can use scripts in bin step by step.

**Stage 0: **
```
alleleCounter -l bed/Select.1000g.snp.bed -b path/to/your/bamfile -o filename.alleleCounter.txt
```
**Stage 1: **
```
perl prepare_ascat_input_after_alleleCounter.pl tumor1.alleleCounter.txt,tumor2.alleleCounter.txt normal1.alleleCounter.txt,normal2.alleleCounter.txt path/of/your/outdir
```

**Stage 2: **
```
Rscript asmultipcf.R path/of/your/outdir(same with stage1)#simplest form to run ASCAT, with matched normal data available, without GC wave correction and all samples female:
or 
Rscript aspcf.R path/of/your/outdir(same with stage1)#run ASCAT with multi-sample segmentation (when shared breakpoints are expected)

```

**Stage 3: **
```
perl add_nummark.pl segment.txt ploidy.txt Tumor_BAF.txt path/of/your/outdir(same with stage1)
```


