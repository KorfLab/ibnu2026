ibnu2026
========

## Notes

- Getting a lightning-specific error with hisat2

## Aligner Test ##

Each aligner runs in its own environment, which is in the `envs` directory.
Install them all as follows:

```
cd envs
ls | xargs -I {} conda env create -f {}
```

At the time of testing, all of the aligners have the same version number in
both macOS-arm64 and linux-64 except bbmap.

| Program      | macOS-arm64 | linux-64   | Notes
|:-------------|:------------|:-----------|:----------------
| bbmap        | 39.77       | 39.79      |
| blast-legacy | 2.2.26      | 2.2.26     |
| bowtie2      | 2.5.5       | 2.5.5      |
| bwa          | 0.7.19      | 0.7.19     |
| gmap         | 2025.07.31  | 2025.07.31 | broken
| hisat2       | 2.2.2       | 2.2.2      |
| minimap2     | 2.30        | 2.30       |
| pblat        | 2.5.1       | 2.5.1      | broken on MacOS
| segemehl     | 0.3.4       | 0.3.4      |
| star         | 2.7.11b     | 2.7.11b    | broken on MacOS
| subread      | 2.1.1       | 2.1.1      |

## Accuracy Experiment ##

- Vary read error from 0-20%
- Calculate %coverage and %missed
- Report resources used

## badread experiment ##

Tried to use badread to generate fake reads but ran into problems

- Negative strand coordinates are wrong (easily fixed from chrom length)
- Read coordinates show one length, but actual length is different
- Hard to say exactly where read is coming from
- Abandoned for now
