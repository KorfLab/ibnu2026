ibnu2026
========

## Aligner Test ##

Each aligner runs in its own environment, which is in the `envs` directory.
Install them all as follows:

```
cd envs
ls | xargs -I {} conda env create -f {}
```

At the time of testing, all of the aligners have the same version number in
both macOS-arm64 and linux-64 except bbmap.

| Program      | macOS-arm64 | linux-64
|:-------------|:------------|:-----------
| bbmap        | 39.77       | 39.79
| blast-legacy | 2.2.26      | 2.2.26
| bowtie2      | 2.5.5       | 2.5.5
| bwa          | 0.7.19      | 0.7.19
| gmap         | 2025.07.31  | 2025.07.31
| hisat2       | 2.2.2       | 2.2.2
| minimap2     | 2.30        | 2.30
| pblat        | 2.5.1       | 2.5.1
| segemehl     | 0.3.4       | 0.3.4
| star         | 2.7.11b     | 2.7.11b
| subread      | 2.1.1       | 2.1.1

