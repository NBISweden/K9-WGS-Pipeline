# K9-WGS-Pipeline
[![Travis status][travis-badge]][travis-link]

Nextflow pipeline for standardised mapping and variant calls on canine genomes

## Usage

```bash
$ nextflow run main.nf
```

## Run on test data

```bash
$ scripts/setup_testdata.sh
$ nextflow run main.nf -profile singularity --fastqDir test-data/test-data-tiny
```

## Authors

- [Johan Viklund](https://github.com/viklund)
- [Jorrit Boekel](https://github.com/glormph)


[travis-badge]: https://api.travis-ci.org/NBISweden/K9-WGS-Pipeline.svg
[travis-link]: https://travis-ci.org/NBISweden/K9-WGS-Pipeline
