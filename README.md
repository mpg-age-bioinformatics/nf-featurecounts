# nf-featurecounts

nextflow differential gene expression pipeline

nextflow homepage: [https://www.nextflow.io](https://www.nextflow.io)

nextflow readthedocs: [https://www.nextflow.io/docs/latest/index.html](https://www.nextflow.io/docs/latest/index.html)

MacOS issues might lead to the need for:
```
export JAVA_HOME=$(/usr/libexec/java_home -v 17)
```

Edit the contents of `project.config` which maps to `nextflow.config` on each respective nextflow (sub)pipe and run:

```
nextflow 1_get_genome
```