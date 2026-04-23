# FragSummarizer

Generate HTML summary reports for [FragPipe](https://fragpipe.nesvilab.org/) results.

FragSummarizer scans a FragPipe output directory and produces a single self-contained HTML report with identification counts, mass error distributions, peptide/charge/missed-cleavage distributions, runtime breakdowns, Percolator feature weights, and MSBooster plots.

## Requirements

- Java 11 or later
- Maven 3.6+ (for building from source)

## Build

```bash
mvn clean package
```

The shaded, runnable JAR is written to `target/FragSummarizer-<version>.jar`.

## Usage

```bash
java -jar target/FragSummarizer-1.0.2.jar -r <fragpipe_results_path>
```

| Flag | Description |
|------|-------------|
| `-r`, `--results_path` | Path to the FragPipe results directory |

The report is written inside the results directory.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
