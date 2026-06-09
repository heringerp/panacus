# Generate an HTML report for your pangenome graph!

> [!TIP]
> You can try this example by downloading this file and running:
> ````bash
> cat html_report.md | sed -n '/^```shell/,/```/p' | sed '/```/d' | bash
> ````

Instead of tab-separated tables, `panacus` supports for many commands also HTML output. The generated report page is interactive and self-contained.

1. Download and unpack the graph:
```shell
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr22.hprc-v1.0-pggb.gfa.gz
gunzip chr22.hprc-v1.0-pggb.gfa.gz
```
2. Prepare file to select subset of paths corresponding to haplotypes (i.e. this tells panacus to remove the references):
```shell
grep '^P' chr22.hprc-v1.0-pggb.gfa | cut -f2 | grep -ve 'grch38\|chm13' > chr22.hprc-v1.0-pggb.paths.haplotypes.txt
```

3. Create the YAML configuration file from which the report should get constructed. This subsets the graph to all non-reference
nodes (see previous step), counts the bps in each node as the "feature", groups paths by sample and tells panacus to run the
histogram (implicitly the histogram will always be calculated, but this tells panacus to show it in the report) and growth curve
analysis.
```shell
echo '- !Gfa
  graph: chr22.hprc-v1.0-pggb.gfa
  subset: chr22.hprc-v1.0-pggb.paths.haplotypes.txt
  count_type: Bp
  grouping: Sample
  analyses:
    - !Hist
    - !Growth
      coverage: 1,2,1,1,1
      quorum: 0,0,1,0.5,0.1' > chr22_pggb.yaml
```

3. Run `panacus report` with 4 threads:
```shell
panacus report -t4 chr22_pggb.yaml > chr22.hprc-v1.0-pggb.histgrowth.html
```

:point_right: :point_right: :point_right: **view the resulting [HTML report here](https://htmlpreview.github.io/?https://github.com/codialab/panacus/blob/main/docs/chr22.hprc-v1.0-pggb.histgrowth.html)!**

> [!TIP]
> If you want to run this for multiple files, e.g. multiple chromosomes, it is very tedious to create multiple YAML files. Instead use keyword replacements to just file in the chromosome number:
> ```shell
> echo '- !Gfa
>   graph: chr{{MY_NUMBER}}.hprc-v1.0-pggb.gfa
>   subset: chr{{MY_NUMBER}}.hprc-v1.0-pggb.paths.haplotypes.txt
>   count_type: Bp
>   grouping: Sample
>   analyses:
>     - !Hist
>     - !Growth
>       coverage: 1,2,1,1,1
>       quorum: 0,0,1,0.5,0.1' > chromosome_wise_pggb.yaml
> ```
> In this case `{{MY_NUMBER}}` will be a placeholder for the number of the chromosome. You can then run this file using:
> ```shell
> panacus report -t4 chromosome_wise_pggb.yaml -r MY_NUMBER=22 > chr22.hprc-v1.0-pggb.histgrowth.html
> ```
> Or in a bash for loop as:
> ```shell
> for i in $(seq 1 22); do panacus report -t4 chromosome_wise_pggb.yaml -r MY_NUMBER="${i}" > chr${i}.hprc-v1.0-pggb.histgrowth.html; done
> ```
> (Remember that you also need to create the subset files for all chromosomes beforehand!)


### Histogram view
![panacus report (coverage histogram) for chr22.hprc-v1.0-pggb.gfa](/docs/chr22.hprc-v1.0-pggb.report.histogram.logscale.highlight.png?raw=true "pangenome report of chr22.hprc-v1.0-pggb.gfa showing coverage histogram in logsacle")


### Growth curve view
![panacus report (pangenome growth) for chr22.hprc-v1.0-pggb.gfa](/docs/chr22.hprc-v1.0-pggb.report.growth.disabled.highlight.png?raw=true "pangenome report of chr22.hprc-v1.0-pggb.gfa showing pangenome growth plots with disabled curves")
