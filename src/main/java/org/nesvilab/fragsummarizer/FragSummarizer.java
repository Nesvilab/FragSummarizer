package org.nesvilab.fragsummarizer;

import tech.tablesaw.api.*;
import tech.tablesaw.columns.Column;
import tech.tablesaw.io.csv.CsvReadOptions;

import java.io.*;
import java.nio.file.*;
import java.time.*;
import java.util.*;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class FragSummarizer {

    private final String resultsPath;

    // Parsed data
    private Table manifestData;
    private Table idNums;
    private Table distributionData;
    private Table ms1MassDf;
    private Table ms2MassDf;
    private LinkedHashMap<String, Double> runtimeDict;
    private LocalDateTime finishTime;
    private Map<String, Table> featuresWeight;
    private List<String> percolatorRawFile;
    private Map<String, Table> singleRunData;
    private Map<String, List<String>> msboosterPlotPaths;

    private Boolean runSpecLib;
    private String latestLogFile;
    private Double ms1Tolerance;
    private Double ms2Tolerance;
    private String ms1Units;
    private String ms2Units;

    // SVG constants
    private static final int ML = 55, MR = 20, MT = 35, MB = 45;
    private static final String GRAY = "#696969";

    public FragSummarizer(String resultsPath) {
        this.resultsPath = resultsPath;
        this.runtimeDict = new LinkedHashMap<>();
        this.featuresWeight = new LinkedHashMap<>();
        this.percolatorRawFile = new ArrayList<>();
        this.singleRunData = new LinkedHashMap<>();
        this.msboosterPlotPaths = new LinkedHashMap<>();
        this.idNums = Table.create("idNums",
                StringColumn.create("Experiment"),
                DoubleColumn.create("PSM"),
                DoubleColumn.create("Peptides"),
                DoubleColumn.create("Proteins"));
        this.distributionData = Table.create("distribution",
                DoubleColumn.create("Peptide Length"),
                DoubleColumn.create("Charge"),
                DoubleColumn.create("Number of Missed Cleavages"),
                StringColumn.create("Exp"));

        readData();
        getRunningTime();
        getMassError();
        getPercolatorFeatures();
        processPsm();
    }

    // ============================================================
    // Tablesaw helpers
    // ============================================================

    /** No-op kept for call-site compatibility. Tablesaw 0.44+ reads CSV columns as StringColumn directly. */
    private static Table convertTextToString(Table t) {
        return t;
    }

    private static String getString(Table t, String colName, int row) {
        return t.column(colName).getString(row);
    }

    private static DoubleColumn asDoubleColumn(Table t, String colName) {
        Column<?> col = t.column(colName);
        if (col instanceof DoubleColumn) return (DoubleColumn) col;
        if (col instanceof IntColumn) return ((IntColumn) col).asDoubleColumn();
        if (col instanceof LongColumn) return ((LongColumn) col).asDoubleColumn();
        if (col instanceof FloatColumn) return ((FloatColumn) col).asDoubleColumn();
        DoubleColumn dc = DoubleColumn.create(colName);
        for (int i = 0; i < col.size(); i++) {
            try { dc.append(Double.parseDouble(col.getString(i))); }
            catch (Exception e) { dc.appendMissing(); }
        }
        return dc;
    }

    // ============================================================
    // Data reading
    // ============================================================

    private void readData() {
        List<String> logFiles = new ArrayList<>();
        File dir = new File(resultsPath);
        String manifestPath = null;

        for (File f : Objects.requireNonNull(dir.listFiles())) {
            if (f.getName().equals("fragpipe-files.fp-manifest")) manifestPath = f.getAbsolutePath();
            if (f.getName().startsWith("log_")) logFiles.add(f.getName());
        }

        Collections.sort(logFiles);
        if (logFiles.isEmpty()) throw new RuntimeException("Log file not found in " + resultsPath);
        latestLogFile = logFiles.get(logFiles.size() - 1);

        if (manifestPath != null) {
            try {
                manifestData = Table.read().csv(CsvReadOptions.builder(manifestPath).separator('\t').header(false).build());
                convertTextToString(manifestData);
                manifestData.column(0).setName("Spectrum File");
                manifestData.column(1).setName("Experiment");
                manifestData.column(2).setName("Bioreplicate");
                manifestData.column(3).setName("Data Type");
                IntColumn runCol = IntColumn.create("Run");
                for (int i = 1; i <= manifestData.rowCount(); i++) runCol.append(i);
                manifestData.addColumns(runCol);
                manifestData = manifestData.reorderColumns("Run", "Spectrum File", "Experiment", "Bioreplicate", "Data Type");
            } catch (Exception e) { throw new RuntimeException("Failed to read manifest", e); }
        }

        try {
            for (String line : Files.readAllLines(Path.of(resultsPath, latestLogFile))) {
                if (line.startsWith("speclibgen.run-speclibgen")) {
                    runSpecLib = !line.split("run-speclibgen=")[1].trim().equals("false");
                    break;
                }
            }
        } catch (IOException e) { throw new RuntimeException("Failed to read log", e); }

        readMsboosterPlots();

        if (runSpecLib != null && runSpecLib && hasDiaData()) {
            readIdsPerRunFromPsm();
        } else {
            appendTable(idNums, readPsm());
            appendTable(idNums, readPeptides());
            appendTable(idNums, readProteins());
        }
    }

    private boolean hasDiaData() {
        if (manifestData == null) return false;
        for (int i = 0; i < manifestData.rowCount(); i++) {
            String val = getString(manifestData, "Data Type", i).trim();
            if (val.equals("1") || val.toUpperCase().contains("DIA")) return true;
        }
        return false;
    }

    private void readMsboosterPlots() {
        Path msboosterDir = Path.of(resultsPath, "MSBooster", "MSBooster_plots");
        if (!Files.exists(msboosterDir)) return;
        try {
            for (File folder : Objects.requireNonNull(msboosterDir.toFile().listFiles())) {
                if (!folder.isDirectory()) continue;
                for (File file : Objects.requireNonNull(folder.listFiles())) {
                    if (!file.getName().endsWith(".png")) continue;
                    String runName = file.getName().contains("edited")
                            ? file.getName().split("_edited")[0]
                            : file.getName().replace(".png", "");
                    boolean shouldInclude = folder.getName().equals("RT_calibration_curves")
                            || file.getName().contains("delta_RT_loess")
                            || file.getName().contains("pred_RT_real_units")
                            || file.getName().contains("unweighted_spectral_entropy");
                    if (shouldInclude) {
                        msboosterPlotPaths.computeIfAbsent(runName, k -> new ArrayList<>())
                                .add(file.getAbsolutePath());
                    }
                }
            }
        } catch (Exception e) { System.err.println("Warning: Failed to read MSBooster plots: " + e.getMessage()); }
    }

    private void readIdsPerRunFromPsm() {
        try {
            Table psmDf = Table.read().csv(CsvReadOptions.builder(Path.of(resultsPath, "psm.tsv").toString()).separator('\t').build());
            convertTextToString(psmDf);

            StringColumn rawFileCol = StringColumn.create("raw_file");
            for (int i = 0; i < psmDf.rowCount(); i++) {
                rawFileCol.append(Path.of(getString(psmDf, "Spectrum", i)).getFileName().toString().split("\\.")[0]);
            }
            psmDf.addColumns(rawFileCol);

            Map<String, String> runToExp = new HashMap<>();
            for (int i = 0; i < manifestData.rowCount(); i++) {
                String basename = Path.of(getString(manifestData, "Spectrum File", i)).getFileName().toString();
                basename = basename.substring(0, basename.lastIndexOf('.'));
                String exp = getString(manifestData, "Experiment", i);
                runToExp.put(basename, (exp == null || exp.trim().isEmpty()) ? basename : exp);
            }

            StringColumn expCol = StringColumn.create("Exp");
            for (int i = 0; i < psmDf.rowCount(); i++)
                expCol.append(runToExp.getOrDefault(rawFileCol.get(i), rawFileCol.get(i)));
            psmDf.addColumns(expCol);

            Table psmCounts = psmDf.countBy(psmDf.stringColumn("Exp"));
            for (int i = 0; i < psmCounts.rowCount(); i++)
                addIdRow(getString(psmCounts, "Exp", i), psmCounts.intColumn("Count").get(i), Double.NaN, Double.NaN);

            if (psmDf.columnNames().contains("Modified Peptide")) {
                Map<String, Set<String>> pepByExp = new HashMap<>();
                for (int i = 0; i < psmDf.rowCount(); i++)
                    pepByExp.computeIfAbsent(expCol.get(i), k -> new HashSet<>()).add(getString(psmDf, "Modified Peptide", i));
                for (int i = 0; i < idNums.rowCount(); i++) {
                    Set<String> s = pepByExp.get(idNums.stringColumn("Experiment").get(i));
                    if (s != null) idNums.doubleColumn("Peptides").set(i, (double) s.size());
                }
            }
            if (psmDf.columnNames().contains("Protein")) {
                Map<String, Set<String>> protByExp = new HashMap<>();
                for (int i = 0; i < psmDf.rowCount(); i++)
                    protByExp.computeIfAbsent(expCol.get(i), k -> new HashSet<>()).add(getString(psmDf, "Protein", i));
                for (int i = 0; i < idNums.rowCount(); i++) {
                    Set<String> s = protByExp.get(idNums.stringColumn("Experiment").get(i));
                    if (s != null) idNums.doubleColumn("Proteins").set(i, (double) s.size());
                }
            }
            if (psmDf.columnNames().contains("Peptide Length") && psmDf.columnNames().contains("Charge")
                    && psmDf.columnNames().contains("Number of Missed Cleavages")) {
                Set<String> seen = new HashSet<>();
                for (int i = 0; i < psmDf.rowCount(); i++) {
                    String key = getString(psmDf, "Spectrum", i) + "\t" + getString(psmDf, "Modified Peptide", i)
                            + "\t" + psmDf.column("Peptide Length").getString(i)
                            + "\t" + psmDf.column("Charge").getString(i)
                            + "\t" + psmDf.column("Number of Missed Cleavages").getString(i);
                    if (seen.add(key)) {
                        distributionData.doubleColumn("Peptide Length").append(toDouble(psmDf, "Peptide Length", i));
                        distributionData.doubleColumn("Charge").append(toDouble(psmDf, "Charge", i));
                        distributionData.doubleColumn("Number of Missed Cleavages").append(toDouble(psmDf, "Number of Missed Cleavages", i));
                        distributionData.stringColumn("Exp").append(expCol.get(i));
                    }
                }
            }
        } catch (Exception e) { System.err.println("Warning: Failed to read PSM per run: " + e.getMessage()); }
    }

    private Table readPsm() {
        Table result = Table.create("psm", StringColumn.create("Experiment"), DoubleColumn.create("PSM"));
        try {
            for (String[] pair : getExpPaths("psm.tsv")) {
                int lines = countLines(pair[1], pair[0]);
                result.stringColumn("Experiment").append(pair[0]);
                result.doubleColumn("PSM").append((double) lines);
            }
        } catch (Exception e) { System.err.println("Warning: Failed to read PSM data: " + e.getMessage()); }
        return result;
    }

    private Table readPeptides() {
        Table result = Table.create("peptides", StringColumn.create("Experiment"), DoubleColumn.create("Peptides"));
        try {
            for (String[] pair : getExpPaths("peptide.tsv")) {
                int lines = countLines(pair[1], "");
                result.stringColumn("Experiment").append(pair[0]);
                result.doubleColumn("Peptides").append((double) lines);
            }
        } catch (Exception e) { System.err.println("Warning: Failed to read peptide data: " + e.getMessage()); }
        return result;
    }

    private Table readProteins() {
        Table result = Table.create("proteins", StringColumn.create("Experiment"), DoubleColumn.create("Proteins"));
        try {
            for (String[] pair : getExpPaths("protein.tsv")) {
                int lines = countLines(pair[1], "");
                result.stringColumn("Experiment").append(pair[0]);
                result.doubleColumn("Proteins").append((double) lines);
            }
        } catch (Exception e) { System.err.println("Warning: Failed to read protein data: " + e.getMessage()); }
        return result;
    }

    /** Returns list of [experimentName, filePath] pairs for the given TSV filename. */
    private List<String[]> getExpPaths(String fileName) {
        List<String[]> result = new ArrayList<>();
        if (isNoExpNoBio() || (runSpecLib != null && runSpecLib)) {
            result.add(new String[]{"One", Path.of(resultsPath, fileName).toString()});
        } else if (isBioNull()) {
            for (String exp : uniqueExperiments())
                result.add(new String[]{exp, Path.of(resultsPath, exp, fileName).toString()});
        } else if (isExpNull()) {
            for (String bio : uniqueBioreplicates())
                result.add(new String[]{bio, Path.of(resultsPath, "exp_" + bio, fileName).toString()});
        } else {
            for (int i = 0; i < manifestData.rowCount(); i++) {
                String key = getString(manifestData, "Experiment", i) + "_" + getString(manifestData, "Bioreplicate", i);
                result.add(new String[]{key, Path.of(resultsPath, key, fileName).toString()});
            }
        }
        return result;
    }

    private void processPsm() {
        try {
            List<String[]> pairs = getExpPaths("psm.tsv");
            for (String[] pair : pairs) {
                if (!Files.exists(Path.of(pair[1]))) continue;
                Table psmDf = Table.read().csv(CsvReadOptions.builder(pair[1]).separator('\t').build());
                convertTextToString(psmDf);

                if (psmDf.columnNames().contains("Retention")) {
                    DoubleColumn retCol = asDoubleColumn(psmDf, "Retention");
                    DoubleColumn retMin = DoubleColumn.create("Retention");
                    for (int i = 0; i < retCol.size(); i++) retMin.append(retCol.get(i) / 60.0);
                    psmDf.replaceColumn("Retention", retMin);
                }

                StringColumn rawFileCol = StringColumn.create("raw_file");
                for (int i = 0; i < psmDf.rowCount(); i++)
                    rawFileCol.append(Path.of(getString(psmDf, "Spectrum", i)).getFileName().toString().split("\\.")[0]);
                psmDf.addColumns(rawFileCol);

                for (String group : new LinkedHashSet<>(rawFileCol.asList()))
                    singleRunData.put(group, psmDf.where(rawFileCol.isEqualTo(group)));
            }
        } catch (Exception e) { System.err.println("Warning: Failed to process PSM data: " + e.getMessage()); }
    }

    private void getRunningTime() {
        try {
            List<String> logLines = Files.readAllLines(Path.of(resultsPath, latestLogFile));
            Map<String, Double> fraggerTime = new LinkedHashMap<>();
            boolean runningTimeRegion = false, mainSearchRegion = false;

            for (String line : logLines) {
                if (!mainSearchRegion) {
                    if (line.startsWith("precursor_true_tolerance = ")) ms1Tolerance = Double.parseDouble(line.split(" = ")[1].trim());
                    if (line.startsWith("fragment_mass_tolerance = ")) ms2Tolerance = Double.parseDouble(line.split(" = ")[1].trim());
                    if (line.startsWith("precursor_true_units = ")) ms1Units = line.split(" = ")[1].trim();
                    if (line.startsWith("fragment_mass_units = ")) ms2Units = line.split(" = ")[1].trim();
                }
                if (line.startsWith("***************************FIRST SEARCH DONE IN"))
                    fraggerTime.put("First Search", Double.parseDouble(line.split("DONE IN ")[1].split(" MIN")[0].trim()));
                if (line.startsWith("************MASS CALIBRATION AND PARAMETER OPTIMIZATION DONE IN")) {
                    fraggerTime.put("Mass Calibration & Param Opt", Double.parseDouble(line.split("DONE IN ")[1].split(" MIN")[0].trim()));
                    mainSearchRegion = true;
                }
                if (line.startsWith("**************************MASS CALIBRATION DONE "))
                    fraggerTime.put("Mass Calibration", Double.parseDouble(line.split("DONE IN ")[1].split(" MIN")[0].trim()));
                if (line.startsWith("***************************MAIN SEARCH DONE IN"))
                    fraggerTime.put("Main Search", Double.parseDouble(line.split("DONE IN ")[1].split(" MIN")[0].trim()));
                if (line.startsWith("Task Runtimes")) { runningTimeRegion = true; continue; }
                if (runningTimeRegion) {
                    String stripped = line.trim();
                    if (stripped.isEmpty()) continue;
                    int lastColon = stripped.lastIndexOf(':');
                    if (lastColon < 0) continue;
                    String task = stripped.substring(0, lastColon);
                    double minutes = Double.parseDouble(stripped.substring(lastColon + 1).trim().split("\\s+")[0]);
                    if (task.equals("MSFragger")) runtimeDict.putAll(fraggerTime);
                    if (minutes != 0 && !task.startsWith("Percolator") && !task.equals("MSFragger")
                            && !task.equals("WorkspaceCleanInit") && !task.equals("WorkspaceClean"))
                        runtimeDict.put(task.trim(), minutes);
                    if (task.startsWith("Percolator")) runtimeDict.merge("Percolator", minutes, Double::sum);
                    if (task.equals("Finalizer Task")) break;
                }
            }
            if ("2".equals(ms1Units) && ms1Tolerance != null) ms1Tolerance *= 1000;
            if ("2".equals(ms2Units) && ms2Tolerance != null) ms2Tolerance *= 1000;
        } catch (Exception e) { System.err.println("Warning: Failed to get running time: " + e.getMessage()); }
    }

    private void getMassError() {
        try {
            List<String> logLines = Files.readAllLines(Path.of(resultsPath, latestLogFile));
            boolean massErrorRegion = false;
            List<String> dataRows = new ArrayList<>();
            percolatorRawFile = new ArrayList<>();
            int massCalPart = 0;

            for (String line : logLines) {
                if (line.startsWith("*********************MASS CALIBRATION AND PARAMETER OPTIMIZATION*******************")
                        || line.startsWith("*********************************MASS CALIBRATION**")) {
                    massErrorRegion = true; massCalPart++;
                }
                if (massErrorRegion && line.trim().matches("^\\d+.*"))
                    dataRows.add(line.trim() + "|" + massCalPart);
                if (line.startsWith("Finding the optimal parameters:") || line.startsWith("**************************MASS CALIBRATION DONE"))
                    massErrorRegion = false;
                if (line.startsWith("Reading tab-delimited input from datafile"))
                    percolatorRawFile.add(line.trim().split("Reading tab-delimited input from datafile ")[1].split("_edited\\.pin")[0]);
            }

            List<Map<String, Object>> ms1Rows = new ArrayList<>(), ms2Rows = new ArrayList<>();
            for (String row : dataRows) {
                String[] parts = row.split("\\|");
                if (parts.length < 6) continue;
                String run = "Part " + parts[5].trim() + ": " + parts[0].trim();
                String[] ms1Old = parts[1].trim().split("\\s+"), ms1New = parts[2].trim().split("\\s+");
                String[] ms2Old = parts[3].trim().split("\\s+"), ms2New = parts[4].trim().split("\\s+");
                if (ms1Old.length >= 2 && ms1New.length >= 2 && ms2Old.length >= 2 && ms2New.length >= 2) {
                    ms1Rows.add(Map.of("Run", run, "Calib", "Old", "Value", Double.parseDouble(ms1Old[0])));
                    ms1Rows.add(Map.of("Run", run, "Calib", "Calibrated", "Value", Double.parseDouble(ms1New[0])));
                    ms2Rows.add(Map.of("Run", run, "Calib", "Old", "Value", Double.parseDouble(ms2Old[0])));
                    ms2Rows.add(Map.of("Run", run, "Calib", "Calibrated", "Value", Double.parseDouble(ms2New[0])));
                }
            }
            ms1MassDf = tableFromMaps("ms1", ms1Rows);
            ms2MassDf = tableFromMaps("ms2", ms2Rows);
        } catch (Exception e) { System.err.println("Warning: Failed to get mass error: " + e.getMessage()); }
    }

    private void getPercolatorFeatures() {
        try {
            String logText = Files.readString(Path.of(resultsPath, latestLogFile));
            Matcher matcher = Pattern.compile(
                    "Learned normalized SVM weights for the 3 cross-validation splits:\\s*\\n(.*?)\\nFound\\s+\\d+",
                    Pattern.DOTALL).matcher(logText);

            int blockIndex = 0;
            while (matcher.find()) {
                String[] lines = matcher.group(1).trim().split("\\n");
                if (lines.length < 2) continue;
                String[] headerCells = lines[0].split("\\t");
                int featureIdx = -1, s1Idx = -1, s2Idx = -1, s3Idx = -1;
                for (int h = 0; h < headerCells.length; h++) {
                    String hdr = headerCells[h].trim();
                    if (hdr.equals("FeatureName")) featureIdx = h;
                    else if (hdr.equals("Split1")) s1Idx = h;
                    else if (hdr.equals("Split2")) s2Idx = h;
                    else if (hdr.equals("Split3")) s3Idx = h;
                }
                if (featureIdx < 0 || s1Idx < 0 || s2Idx < 0 || s3Idx < 0) continue;

                List<String> names = new ArrayList<>();
                List<Double> means = new ArrayList<>();
                for (int i = 1; i < lines.length; i++) {
                    if (lines[i].trim().isEmpty()) continue;
                    String[] cells = lines[i].split("\\t");
                    int maxIdx = Math.max(featureIdx, Math.max(s1Idx, Math.max(s2Idx, s3Idx)));
                    if (cells.length <= maxIdx) continue;
                    String name = cells[featureIdx].trim();
                    if (name.contains("m0")) continue;
                    double s1 = Double.parseDouble(cells[s1Idx].trim());
                    double s2 = Double.parseDouble(cells[s2Idx].trim());
                    double s3 = Double.parseDouble(cells[s3Idx].trim());
                    if (s1 == 0 && s2 == 0 && s3 == 0) continue;
                    names.add(name);
                    means.add(Math.round(((s1 + s2 + s3) / 3.0) * 10000.0) / 10000.0);
                }

                Integer[] idx = new Integer[names.size()];
                for (int i = 0; i < idx.length; i++) idx[i] = i;
                Arrays.sort(idx, (a, b) -> Double.compare(Math.abs(means.get(b)), Math.abs(means.get(a))));

                Table features = Table.create("features", StringColumn.create("FeatureName"), DoubleColumn.create("Mean Weight"));
                for (int ii : idx) {
                    features.stringColumn("FeatureName").append(names.get(ii));
                    features.doubleColumn("Mean Weight").append(means.get(ii));
                }
                if (blockIndex < percolatorRawFile.size())
                    featuresWeight.put(percolatorRawFile.get(blockIndex), features);
                blockIndex++;
            }

            // Parse finish time
            List<Duration> durations = new ArrayList<>();
            Matcher m1 = Pattern.compile("time=\"(\\d{2}):(\\d{2}):(\\d{2})\"").matcher(logText);
            while (m1.find()) durations.add(Duration.ofHours(Integer.parseInt(m1.group(1))).plusMinutes(Integer.parseInt(m1.group(2))).plusSeconds(Integer.parseInt(m1.group(3))));
            Matcher m2 = Pattern.compile("\\[\\d{4}/\\d{2}/\\d{2} (\\d{2}):(\\d{2}):(\\d{2})\\]").matcher(logText);
            while (m2.find()) durations.add(Duration.ofHours(Integer.parseInt(m2.group(1))).plusMinutes(Integer.parseInt(m2.group(2))).plusSeconds(Integer.parseInt(m2.group(3))));

            if (!durations.isEmpty()) {
                LocalDate creationDate = Files.getLastModifiedTime(Path.of(resultsPath, latestLogFile)).toInstant().atZone(ZoneId.systemDefault()).toLocalDate();
                finishTime = LocalDateTime.of(creationDate, LocalTime.MIDNIGHT).plus(durations.stream().max(Comparator.naturalOrder()).orElse(Duration.ZERO));
            }
        } catch (Exception e) { System.err.println("Warning: Failed to parse percolator features: " + e.getMessage()); }
    }

    // ============================================================
    // HTML report generation (static SVG, no JavaScript)
    // ============================================================

    public void generateReport() {
        StringBuilder html = new StringBuilder();
        html.append(htmlHeader());
        html.append(overviewTableHtml());
        html.append(ganttChartHtml());
        html.append(overallStatsHtml());

        for (String runName : singleRunData.keySet())
            html.append(runPageHtml(runName));

        for (String runName : singleRunData.keySet())
            if (msboosterPlotPaths.containsKey(runName))
                html.append(msboosterPageHtml(runName));

        html.append("</body>\n</html>\n");

        try {
            String outputPath = Path.of(resultsPath, "fragpipe-report.html").toString();
            Files.writeString(Path.of(outputPath), html.toString());
            System.out.println("Report generated successfully: " + outputPath);
        } catch (IOException e) {
            System.err.println("Failed to write report: " + e.getMessage());
        }
    }

    private String htmlHeader() {
        return "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n"
                + "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
                + "<title>FragPipe Report</title>\n"
                + "<style>\n"
                + "* { box-sizing: border-box; margin: 0; padding: 0; }\n"
                + "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;\n"
                + "  background: #f5f5f5; color: #333; }\n"
                + ".page { background: #fff; max-width: 1200px; margin: 24px auto; padding: 32px;\n"
                + "  border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }\n"
                + ".page h2 { font-size: 20px; margin-bottom: 16px; color: #1a1a1a;\n"
                + "  border-bottom: 2px solid #A0E7E5; padding-bottom: 8px; }\n"
                + "table.overview { width: 100%; border-collapse: collapse; font-size: 13px; }\n"
                + "table.overview th { background: #A0E7E5; padding: 8px 10px; text-align: left;\n"
                + "  font-weight: 600; border-bottom: 2px solid #fff; }\n"
                + "table.overview td { padding: 7px 10px; border-bottom: 1px solid #e8e8e8; }\n"
                + "table.overview tr:hover td { background: #f8fffe; }\n"
                + ".grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }\n"
                + ".grid3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 8px; }\n"
                + ".chart-cell svg { width: 100%; height: auto; display: block; }\n"
                + ".total-badge { display: inline-block; background: #fff3f3; border: 1px solid #e74c3c;\n"
                + "  border-radius: 6px; padding: 4px 12px; font-weight: 600; color: #e74c3c; float: right; }\n"
                + ".img-row { display: flex; flex-wrap: wrap; gap: 12px; justify-content: center; align-items: flex-start; }\n"
                + ".img-row img { max-width: 48%; height: auto; object-fit: contain; border: 1px solid #eee; border-radius: 4px; }\n"
                + "@media print { .page { box-shadow: none; margin: 0; page-break-after: always; } }\n"
                + "</style>\n</head>\n<body>\n";
    }

    // ---- Overview table ----

    private String overviewTableHtml() {
        if (manifestData == null) return "";
        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"page\">\n<h2>Experiment Design</h2>\n");
        sb.append("<table class=\"overview\">\n<tr>");
        for (String col : manifestData.columnNames()) sb.append("<th>").append(esc(col)).append("</th>");
        sb.append("</tr>\n");
        for (int r = 0; r < manifestData.rowCount(); r++) {
            sb.append("<tr>");
            for (int c = 0; c < manifestData.columnCount(); c++)
                sb.append("<td>").append(esc(manifestData.column(c).getString(r))).append("</td>");
            sb.append("</tr>\n");
        }
        sb.append("</table>\n</div>\n");
        return sb.toString();
    }

    // ---- Gantt chart ----

    private String ganttChartHtml() {
        if (runtimeDict == null || runtimeDict.isEmpty()) return "";
        double total = runtimeDict.values().stream().mapToDouble(Double::doubleValue).sum();

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"page\">\n<h2>Workflow Processing Time")
                .append("<span class=\"total-badge\">Total: ").append(String.format("%.2f", total)).append(" min</span></h2>\n");

        if (finishTime != null) {
            sb.append(svgGantt());
        } else {
            // Fallback: simple duration bar chart when no finish time available
            List<String> tasks = new ArrayList<>(runtimeDict.keySet());
            Collections.reverse(tasks);
            List<Double> values = tasks.stream().map(runtimeDict::get).collect(Collectors.toList());
            sb.append(svgHBar(tasks, values, "", "Duration (min)", "#6B8EDB", 700));
        }
        sb.append("</div>\n");
        return sb.toString();
    }

    /** Proper Gantt chart: tasks positioned on a time axis (x = clock time HH:MM) */
    private String svgGantt() {
        // Compute start/end times working backwards from finishTime
        List<String> taskNames = new ArrayList<>();
        List<LocalDateTime> starts = new ArrayList<>();
        List<LocalDateTime> ends = new ArrayList<>();
        List<Double> durations = new ArrayList<>();

        LocalDateTime currentEnd = finishTime;
        List<Map.Entry<String, Double>> entries = new ArrayList<>(runtimeDict.entrySet());
        Collections.reverse(entries);
        for (Map.Entry<String, Double> e : entries) {
            double dur = e.getValue();
            LocalDateTime taskStart = currentEnd.minusSeconds((long) (dur * 60));
            taskNames.add(0, e.getKey());
            starts.add(0, taskStart);
            ends.add(0, currentEnd);
            durations.add(0, dur);
            currentEnd = taskStart;
        }

        int n = taskNames.size();
        LocalDateTime globalStart = starts.get(0);
        LocalDateTime globalEnd = ends.get(n - 1);
        long totalSeconds = java.time.Duration.between(globalStart, globalEnd).getSeconds();
        if (totalSeconds <= 0) totalSeconds = 1;

        // SVG dimensions
        int svgW = 800;
        int labelW = 200;
        int pw = svgW - labelW - MR;
        int hBarH = 22, hGap = 6;
        int svgH = MT + n * (hBarH + hGap) + MB + 15;
        int ph = n * (hBarH + hGap);

        // Time tick interval (minutes)
        double totalMin = totalSeconds / 60.0;
        int interval;
        if (totalMin <= 10) interval = 1;
        else if (totalMin <= 30) interval = 5;
        else if (totalMin <= 120) interval = 15;
        else if (totalMin <= 360) interval = 30;
        else interval = 60;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" style=\"width:100%%;max-width:%dpx\">\n", svgW, svgH, svgW));

        // Vertical grid lines at time intervals
        // Round globalStart down to the nearest interval
        LocalDateTime tickTime = globalStart.withSecond(0).withNano(0);
        int startMin = tickTime.getHour() * 60 + tickTime.getMinute();
        startMin = (startMin / interval) * interval;
        tickTime = tickTime.withHour(startMin / 60).withMinute(startMin % 60);
        if (tickTime.isBefore(globalStart)) tickTime = tickTime.plusMinutes(interval);

        while (!tickTime.isAfter(globalEnd)) {
            long secFromStart = java.time.Duration.between(globalStart, tickTime).getSeconds();
            double x = labelW + (double) secFromStart / totalSeconds * pw;
            sb.append(String.format("<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e0e0e0\" stroke-width=\"1\"/>\n", x, MT, x, MT + ph));
            sb.append(svgText(x, MT + ph + 15, String.format("%02d:%02d", tickTime.getHour(), tickTime.getMinute()), "middle", 9, false));
            tickTime = tickTime.plusMinutes(interval);
        }

        // Task bars (top-to-bottom = chronological)
        for (int i = 0; i < n; i++) {
            long sStart = java.time.Duration.between(globalStart, starts.get(i)).getSeconds();
            long sEnd = java.time.Duration.between(globalStart, ends.get(i)).getSeconds();
            double x1 = labelW + (double) sStart / totalSeconds * pw;
            double x2 = labelW + (double) sEnd / totalSeconds * pw;
            double y = MT + i * (hBarH + hGap);
            double barW = Math.max(x2 - x1, 1);

            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%d\" fill=\"#6B8EDB\" rx=\"2\"/>\n",
                    x1, y, barW, hBarH));

            // Duration label on bar
            String durLabel = String.format("%.2f min", durations.get(i));
            double textX = x1 + barW / 2;
            if (barW > 80) {
                // Label inside bar
                sb.append(String.format("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"8\" font-weight=\"bold\" fill=\"#000\">%s</text>\n",
                        textX, y + hBarH / 2.0 + 3, esc(durLabel)));
            } else {
                // Label outside bar (right)
                sb.append(String.format("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"start\" font-size=\"8\" font-weight=\"bold\">%s</text>\n",
                        x2 + 3, y + hBarH / 2.0 + 3, esc(durLabel)));
            }

            // Task name (left)
            sb.append(svgText(labelW - 5, y + hBarH / 2.0 + 4, taskNames.get(i), "end", 10, false));
        }

        // X axis line
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", labelW, MT + ph, labelW + pw, MT + ph));

        // X axis label
        sb.append(svgText(labelW + pw / 2.0, svgH - 3, "Time", "middle", 11, false));

        sb.append("</svg>\n");
        return sb.toString();
    }

    // ---- Overall statistics page ----

    private String overallStatsHtml() {
        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"page\">\n<h2>Overall Statistics</h2>\n");

        // Row 1: Charge + Peptide Length
        sb.append("<div class=\"grid2\">\n");
        sb.append(svgCategoryBar(asDoubleColumn(distributionData, "Charge"), "Charge", "Charge States", "# PSM", 400, 300));
        sb.append(svgHistogram(toDoubleArrayNoNaN(asDoubleColumn(distributionData, "Peptide Length")), "Peptide Length", "Peptide Length", "# PSM", 400, 300));
        sb.append("</div>\n");

        // Row 2: PSM, Peptides, Proteins
        // Group idNums by Experiment, taking the first non-NaN value per metric
        // (matches Python's groupby('Experiment').agg('first'))
        Table grouped = groupIdNums();
        if (grouped.rowCount() > 0) {
            sb.append("<div class=\"grid3\">\n");
            for (String metric : new String[]{"PSM", "Peptides", "Proteins"}) {
                if (grouped.columnNames().contains(metric))
                    sb.append(svgExperimentBar(grouped.stringColumn("Experiment"), asDoubleColumn(grouped, metric),
                            metric + " IDs", "Experiments", "# " + metric, 400, 300));
            }
            sb.append("</div>\n");
        }

        // Row 3: Missed Cleavage, MS1 error, MS2 error
        sb.append("<div class=\"grid3\">\n");
        sb.append(svgCategoryBar(asDoubleColumn(distributionData, "Number of Missed Cleavages"), "Missed Cleavage", "Missed Cleavages", "# PSM", 400, 300));
        sb.append(svgMassError(ms1MassDf, "MS1 Mass Error", 400, 300));
        sb.append(svgMassError(ms2MassDf, "MS2 Mass Error", 400, 300));
        sb.append("</div>\n");

        sb.append("</div>\n");
        return sb.toString();
    }

    // ---- Per-run page ----

    private String runPageHtml(String runName) {
        if (!singleRunData.containsKey(runName)) return "";
        Table data = singleRunData.get(runName);

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"page\">\n<h2>").append(esc(runName)).append(" Statistics</h2>\n");

        // Row 1
        sb.append("<div class=\"grid3\">\n");
        if (hasCol(data, "Peptide Length")) sb.append(svgHistogram(toDoubleArrayNoNaN(asDoubleColumn(data, "Peptide Length")), "Peptide Length Distribution", "Peptide Length", "# PSM", 400, 300));
        if (hasCol(data, "Number of Missed Cleavages")) sb.append(svgCategoryBar(asDoubleColumn(data, "Number of Missed Cleavages"), "Missed Cleavages", "Missed Cleavages", "# PSM", 400, 300));
        if (hasCol(data, "Charge")) sb.append(svgCategoryBar(asDoubleColumn(data, "Charge"), "Charge", "Charge", "# PSM", 400, 300));
        sb.append("</div>\n");

        // Row 2
        sb.append("<div class=\"grid3\">\n");
        if (hasCol(data, "Calculated M/Z")) sb.append(svgHistogram(toDoubleArrayNoNaN(asDoubleColumn(data, "Calculated M/Z")), "M/Z Distribution", "M/Z", "# PSM", 400, 300));
        if (hasCol(data, "Retention")) sb.append(svgHistogram(toDoubleArrayNoNaN(asDoubleColumn(data, "Retention")), "Retention Time Distribution", "Retention Time (min)", "# PSM", 400, 300));
        if (hasCol(data, "Hyperscore")) sb.append(svgHistogram(toDoubleArrayNoNaN(asDoubleColumn(data, "Hyperscore")), "Hyperscore Distribution", "Hyperscore", "# PSM", 400, 300));
        sb.append("</div>\n");

        // Row 3
        sb.append("<div class=\"grid3\">\n");
        if (hasCol(data, "Delta Mass")) sb.append(svgHistogram(toDoubleArrayNoNaN(asDoubleColumn(data, "Delta Mass")), "Delta Mass Dis", "Delta Mass (Da)", "# PSM", 400, 300));
        if (hasCol(data, "Calculated M/Z") && hasCol(data, "Delta Mass"))
            sb.append(svgScatter(asDoubleColumn(data, "Calculated M/Z"), asDoubleColumn(data, "Delta Mass"), 0.5, "M/Z vs Delta Mass", "M/Z", "Delta Mass (Da)", 400, 300));
        if (hasCol(data, "Retention") && hasCol(data, "Calculated M/Z"))
            sb.append(svgScatter(asDoubleColumn(data, "Retention"), asDoubleColumn(data, "Calculated M/Z"), -1, "RT vs M/Z", "Retention Time (min)", "M/Z", 400, 300));
        sb.append("</div>\n");

        sb.append("</div>\n");
        return sb.toString();
    }

    // ---- MSBooster + Percolator features page ----

    private String msboosterPageHtml(String runName) {
        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"page\">\n<h2>MSBooster &amp; Percolator Features &mdash; ").append(esc(runName)).append("</h2>\n");

        // Percolator features horizontal bar
        if (featuresWeight.containsKey(runName)) {
            Table df = featuresWeight.get(runName);
            if (df.rowCount() > 0) {
                List<String> names = new ArrayList<>();
                List<Double> weights = new ArrayList<>();
                List<String> colors = new ArrayList<>();
                for (int i = df.rowCount() - 1; i >= 0; i--) {
                    names.add(df.stringColumn("FeatureName").get(i));
                    double w = df.doubleColumn("Mean Weight").get(i);
                    weights.add(w);
                    colors.add(w > 0 ? "#2ecc71" : "#e74c3c");
                }
                sb.append(svgHBarColored(names, weights, colors, "Percolator Feature Weights", "Mean Weight", 700));
            }
        }

        // MSBooster images
        List<String> paths = msboosterPlotPaths.get(runName);
        if (paths != null && !paths.isEmpty()) {
            sb.append("<div class=\"img-row\">\n");
            for (String imgPath : paths) {
                String b64 = imageToBase64(imgPath);
                if (b64 != null) sb.append("<img src=\"data:image/png;base64,").append(b64).append("\" alt=\"MSBooster plot\">\n");
            }
            sb.append("</div>\n");
        }

        sb.append("</div>\n");
        return sb.toString();
    }

    // ============================================================
    // SVG chart builders
    // ============================================================

    private String svgHistogram(double[] vals, String title, String xLabel, String yLabel, int w, int h) {
        if (vals.length == 0) return "<div class=\"chart-cell\"></div>\n";
        int pw = w - ML - MR, ph = h - MT - MB;

        double dataMin = Double.MAX_VALUE, dataMax = -Double.MAX_VALUE;
        for (double v : vals) { if (v < dataMin) dataMin = v; if (v > dataMax) dataMax = v; }
        if (dataMin == dataMax) { dataMin -= 0.5; dataMax += 0.5; }

        // Detect integer-valued data: if every value is a whole number, use 1 bin per integer
        boolean isInteger = true;
        for (double v : vals) { if (v != Math.floor(v)) { isInteger = false; break; } }

        int nbins;
        double binW;
        double binMin; // left edge of the first bin
        if (isInteger && (dataMax - dataMin) < 200) {
            // One bin per integer value — no gaps
            int lo = (int) dataMin, hi = (int) dataMax;
            nbins = hi - lo + 1;
            binW = 1.0;
            binMin = lo;
        } else {
            nbins = Math.min(50, Math.max(10, vals.length / 20));
            binW = (dataMax - dataMin) / nbins;
            binMin = dataMin;
        }

        int[] counts = new int[nbins];
        for (double v : vals) {
            int b = (int) ((v - binMin) / binW);
            if (b >= nbins) b = nbins - 1;
            if (b < 0) b = 0;
            counts[b]++;
        }
        int maxCount = 1;
        for (int c : counts) if (c > maxCount) maxCount = c;

        double[] ys = niceScale(0, maxCount, 5);
        double yMax = ys[1], yStep = ys[2];
        if (yMax == 0) yMax = 1;

        double totalRange = nbins * binW;

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"chart-cell\">\n");
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", w, h));

        // Title
        sb.append(svgText(w / 2.0, 20, title, "middle", 13, true));

        // Grid lines
        for (double tick = yStep; tick <= yMax; tick += yStep) {
            double y = MT + ph - (tick / yMax * ph);
            sb.append(String.format("<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n", ML, y, ML + pw, y));
        }

        // Bars
        double barPx = (double) pw / nbins;
        for (int i = 0; i < nbins; i++) {
            if (counts[i] == 0) continue;
            double x = ML + i * barPx;
            double barH = counts[i] / yMax * ph;
            double y = MT + ph - barH;
            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n",
                    x, y, barPx, barH, GRAY));
        }

        // Axes
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT, ML, MT + ph));
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT + ph, ML + pw, MT + ph));

        // Y axis ticks
        for (double tick = 0; tick <= yMax + yStep * 0.01; tick += yStep) {
            double y = MT + ph - (tick / yMax * ph);
            sb.append(svgText(ML - 5, y + 4, fmtTick(tick), "end", 9, false));
        }

        // X axis ticks
        double[] xs = niceScale(binMin, binMin + totalRange, 5);
        for (double tick = xs[0]; tick <= xs[1] + xs[2] * 0.01; tick += xs[2]) {
            if (tick < binMin || tick > binMin + totalRange) continue;
            double x = ML + (tick - binMin) / totalRange * pw;
            sb.append(svgText(x, MT + ph + 15, fmtTick(tick), "middle", 9, false));
        }

        // Labels
        sb.append(svgText(ML + pw / 2.0, h - 5, xLabel, "middle", 11, false));
        sb.append(svgRotText(12, MT + ph / 2.0, yLabel, 11));

        sb.append("</svg>\n</div>\n");
        return sb.toString();
    }

    private String svgCategoryBar(DoubleColumn data, String title, String xLabel, String yLabel, int w, int h) {
        double[] vals = toDoubleArrayNoNaN(data);
        if (vals.length == 0) return "<div class=\"chart-cell\"></div>\n";
        int pw = w - ML - MR, ph = h - MT - MB;

        Map<Integer, Integer> counts = new TreeMap<>();
        for (double v : vals) counts.merge((int) v, 1, Integer::sum);

        List<Integer> keys = new ArrayList<>(counts.keySet());
        int maxCount = counts.values().stream().mapToInt(Integer::intValue).max().orElse(1);
        double[] ys = niceScale(0, maxCount, 5);
        double yMax = ys[1], yStep = ys[2];
        if (yMax == 0) yMax = 1;

        int n = keys.size();
        double groupW = (double) pw / n;
        double barW = groupW * 0.7;

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"chart-cell\">\n");
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", w, h));
        sb.append(svgText(w / 2.0, 20, title, "middle", 13, true));

        // Grid
        for (double tick = yStep; tick <= yMax; tick += yStep) {
            double y = MT + ph - (tick / yMax * ph);
            sb.append(String.format("<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n", ML, y, ML + pw, y));
        }

        // Bars + labels
        for (int i = 0; i < n; i++) {
            int count = counts.get(keys.get(i));
            double x = ML + i * groupW + (groupW - barW) / 2;
            double barH = count / yMax * ph;
            double y = MT + ph - barH;
            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n", x, y, barW, barH, GRAY));
            sb.append(svgText(ML + i * groupW + groupW / 2, y - 4, String.valueOf(count), "middle", 9, false));
            sb.append(svgText(ML + i * groupW + groupW / 2, MT + ph + 15, String.valueOf(keys.get(i)), "middle", 9, false));
        }

        // Axes
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT, ML, MT + ph));
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT + ph, ML + pw, MT + ph));

        // Y ticks
        for (double tick = 0; tick <= yMax + yStep * 0.01; tick += yStep) {
            double y = MT + ph - (tick / yMax * ph);
            sb.append(svgText(ML - 5, y + 4, fmtTick(tick), "end", 9, false));
        }

        // Labels
        sb.append(svgText(ML + pw / 2.0, h - 5, xLabel, "middle", 11, false));
        sb.append(svgRotText(12, MT + ph / 2.0, yLabel, 11));

        sb.append("</svg>\n</div>\n");
        return sb.toString();
    }

    private String svgExperimentBar(StringColumn exps, DoubleColumn vals, String title, String xLabel, String yLabel, int w, int h) {
        if (exps.size() == 0) return "<div class=\"chart-cell\"></div>\n";
        // Use more bottom margin for rotated labels
        int mbLocal = 70;
        int pw = w - ML - MR, ph = h - MT - mbLocal;

        List<String> names = new ArrayList<>();
        List<Double> values = new ArrayList<>();
        for (int i = 0; i < exps.size(); i++) {
            names.add(exps.get(i));
            values.add(vals.isMissing(i) ? 0.0 : vals.get(i));
        }

        double maxVal = values.stream().mapToDouble(Double::doubleValue).max().orElse(1);
        double[] ys = niceScale(0, maxVal, 5);
        double yMax = ys[1], yStep = ys[2];
        if (yMax == 0) yMax = 1;

        int n = names.size();
        double groupW = (double) pw / n;
        double barW = groupW * 0.7;

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"chart-cell\">\n");
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", w, h));
        sb.append(svgText(w / 2.0, 20, title, "middle", 13, true));

        // Grid
        for (double tick = yStep; tick <= yMax; tick += yStep) {
            double y = MT + ph - (tick / yMax * ph);
            sb.append(String.format("<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n", ML, y, ML + pw, y));
        }

        // Bars
        for (int i = 0; i < n; i++) {
            double v = values.get(i);
            double x = ML + i * groupW + (groupW - barW) / 2;
            double barH = v / yMax * ph;
            double y = MT + ph - barH;
            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"%s\"/>\n", x, y, barW, barH, GRAY));
            sb.append(svgText(ML + i * groupW + groupW / 2, y - 4, String.valueOf((int) Math.round(v)), "middle", 9, true));

            // Rotated x label
            double lx = ML + i * groupW + groupW / 2;
            double ly = MT + ph + 8;
            String label = names.get(i).length() > 20 ? names.get(i).substring(0, 20) + ".." : names.get(i);
            sb.append(String.format("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-size=\"8\" transform=\"rotate(-40,%.1f,%.1f)\">%s</text>\n",
                    lx, ly, lx, ly, esc(label)));
        }

        // Axes
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT, ML, MT + ph));
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT + ph, ML + pw, MT + ph));

        // Y ticks
        for (double tick = 0; tick <= yMax + yStep * 0.01; tick += yStep) {
            double y = MT + ph - (tick / yMax * ph);
            sb.append(svgText(ML - 5, y + 4, fmtTick(tick), "end", 9, false));
        }

        sb.append(svgRotText(12, MT + ph / 2.0, yLabel, 11));

        sb.append("</svg>\n</div>\n");
        return sb.toString();
    }

    private String svgScatter(DoubleColumn xCol, DoubleColumn yCol, double yThreshold, String title, String xLabel, String yLabel, int w, int h) {
        List<Double> xList = new ArrayList<>(), yList = new ArrayList<>();
        for (int i = 0; i < xCol.size(); i++) {
            if (xCol.isMissing(i) || yCol.isMissing(i)) continue;
            if (yThreshold > 0 && Math.abs(yCol.get(i)) >= yThreshold) continue;
            xList.add(xCol.get(i));
            yList.add(yCol.get(i));
        }
        if (xList.isEmpty()) return "<div class=\"chart-cell\"></div>\n";

        // Downsample if too many points
        if (xList.size() > 5000) {
            int step = xList.size() / 5000;
            List<Double> sx = new ArrayList<>(), sy = new ArrayList<>();
            for (int i = 0; i < xList.size(); i += step) { sx.add(xList.get(i)); sy.add(yList.get(i)); }
            xList = sx; yList = sy;
        }

        int pw = w - ML - MR, ph = h - MT - MB;
        double xMin = xList.stream().mapToDouble(Double::doubleValue).min().orElse(0);
        double xMax = xList.stream().mapToDouble(Double::doubleValue).max().orElse(1);
        double yMin = yList.stream().mapToDouble(Double::doubleValue).min().orElse(0);
        double yMax = yList.stream().mapToDouble(Double::doubleValue).max().orElse(1);
        if (xMin == xMax) { xMin -= 0.5; xMax += 0.5; }
        if (yMin == yMax) { yMin -= 0.5; yMax += 0.5; }

        double[] xs = niceScale(xMin, xMax, 5);
        double[] ys = niceScale(yMin, yMax, 5);

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"chart-cell\">\n");
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", w, h));
        sb.append(svgText(w / 2.0, 20, title, "middle", 13, true));

        // Grid
        for (double tick = ys[0]; tick <= ys[1] + ys[2] * 0.01; tick += ys[2]) {
            double y = MT + ph - ((tick - ys[0]) / (ys[1] - ys[0]) * ph);
            sb.append(String.format("<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n", ML, y, ML + pw, y));
        }

        // Points
        sb.append("<g fill=\"#808080\" fill-opacity=\"0.4\">\n");
        double xRange = xMax - xMin, yRange = yMax - yMin;
        for (int i = 0; i < xList.size(); i++) {
            double px = ML + (xList.get(i) - xMin) / xRange * pw;
            double py = MT + ph - (yList.get(i) - yMin) / yRange * ph;
            sb.append(String.format("<circle cx=\"%.1f\" cy=\"%.1f\" r=\"1.5\"/>\n", px, py));
        }
        sb.append("</g>\n");

        // Axes
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT, ML, MT + ph));
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT + ph, ML + pw, MT + ph));

        // Y ticks
        for (double tick = ys[0]; tick <= ys[1] + ys[2] * 0.01; tick += ys[2]) {
            double y = MT + ph - ((tick - ys[0]) / (ys[1] - ys[0]) * ph);
            sb.append(svgText(ML - 5, y + 4, fmtTick(tick), "end", 9, false));
        }

        // X ticks
        for (double tick = xs[0]; tick <= xs[1] + xs[2] * 0.01; tick += xs[2]) {
            if (tick < xMin - xs[2] * 0.5 || tick > xMax + xs[2] * 0.5) continue;
            double x = ML + (tick - xMin) / xRange * pw;
            sb.append(svgText(x, MT + ph + 15, fmtTick(tick), "middle", 9, false));
        }

        sb.append(svgText(ML + pw / 2.0, h - 5, xLabel, "middle", 11, false));
        sb.append(svgRotText(12, MT + ph / 2.0, yLabel, 11));

        sb.append("</svg>\n</div>\n");
        return sb.toString();
    }

    private String svgMassError(Table massDf, String title, int w, int h) {
        if (massDf == null || massDf.rowCount() == 0) return "<div class=\"chart-cell\"></div>\n";
        int mbLocal = 70;
        int pw = w - ML - MR, ph = h - MT - mbLocal;

        List<String> runs = massDf.stringColumn("Run").unique().asList().stream().sorted().collect(Collectors.toList());
        List<Double> oldVals = new ArrayList<>(), calVals = new ArrayList<>();
        for (String run : runs) {
            Table rd = massDf.where(massDf.stringColumn("Run").isEqualTo(run));
            double ov = 0, cv = 0;
            for (int i = 0; i < rd.rowCount(); i++) {
                if ("Old".equals(rd.stringColumn("Calib").get(i))) ov = rd.doubleColumn("Value").get(i);
                else cv = rd.doubleColumn("Value").get(i);
            }
            oldVals.add(ov); calVals.add(cv);
        }

        double allMin = Double.MAX_VALUE, allMax = -Double.MAX_VALUE;
        for (double v : oldVals) { allMin = Math.min(allMin, v); allMax = Math.max(allMax, v); }
        for (double v : calVals) { allMin = Math.min(allMin, v); allMax = Math.max(allMax, v); }
        double[] ys = niceScale(Math.min(allMin, 0), Math.max(allMax, 0), 5);
        double yLo = ys[0], yHi = ys[1], yStep = ys[2];
        double yRange = yHi - yLo;
        if (yRange == 0) yRange = 1;

        int n = runs.size();
        double groupW = (double) pw / n;
        double barW = groupW * 0.35;

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"chart-cell\">\n");
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\">\n", w, h));
        sb.append(svgText(w / 2.0, 20, title, "middle", 13, true));

        // Grid
        for (double tick = yLo; tick <= yHi + yStep * 0.01; tick += yStep) {
            double y = MT + ph - ((tick - yLo) / yRange * ph);
            sb.append(String.format("<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n", ML, y, ML + pw, y));
        }

        // Zero line
        if (yLo < 0 && yHi > 0) {
            double y0 = MT + ph - ((0 - yLo) / yRange * ph);
            sb.append(String.format("<line x1=\"%d\" y1=\"%.1f\" x2=\"%d\" y2=\"%.1f\" stroke=\"#999\" stroke-width=\"1\" stroke-dasharray=\"4,2\"/>\n", ML, y0, ML + pw, y0));
        }

        double baseline = MT + ph - ((0 - yLo) / yRange * ph);
        for (int i = 0; i < n; i++) {
            double cx = ML + i * groupW + groupW / 2;
            // Old bar
            double ov = oldVals.get(i);
            double oy = MT + ph - ((ov - yLo) / yRange * ph);
            double oh = Math.abs(oy - baseline);
            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#808080\"/>\n",
                    cx - barW - 1, Math.min(oy, baseline), barW, oh));
            // Calibrated bar
            double cv = calVals.get(i);
            double cy = MT + ph - ((cv - yLo) / yRange * ph);
            double ch = Math.abs(cy - baseline);
            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"#1a1a1a\"/>\n",
                    cx + 1, Math.min(cy, baseline), barW, ch));

            // Run label (rotated)
            double lx = cx;
            double ly = MT + ph + 8;
            String label = runs.get(i).length() > 15 ? runs.get(i).substring(0, 15) + ".." : runs.get(i);
            sb.append(String.format("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" font-size=\"7\" transform=\"rotate(-40,%.1f,%.1f)\">%s</text>\n",
                    lx, ly, lx, ly, esc(label)));
        }

        // Axes
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT, ML, MT + ph));
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", ML, MT + ph, ML + pw, MT + ph));

        // Y ticks
        for (double tick = yLo; tick <= yHi + yStep * 0.01; tick += yStep) {
            double y = MT + ph - ((tick - yLo) / yRange * ph);
            sb.append(svgText(ML - 5, y + 4, fmtTick(tick), "end", 9, false));
        }

        // Legend
        sb.append(String.format("<rect x=\"%d\" y=\"%d\" width=\"10\" height=\"10\" fill=\"#808080\"/>\n", ML + pw - 120, MT + 2));
        sb.append(svgText(ML + pw - 107, MT + 11, "Before", "start", 9, false));
        sb.append(String.format("<rect x=\"%d\" y=\"%d\" width=\"10\" height=\"10\" fill=\"#1a1a1a\"/>\n", ML + pw - 60, MT + 2));
        sb.append(svgText(ML + pw - 47, MT + 11, "Calibrated", "start", 9, false));

        sb.append(svgRotText(12, MT + ph / 2.0, "Median Mass Error", 11));

        sb.append("</svg>\n</div>\n");
        return sb.toString();
    }

    /** Horizontal bar chart with uniform color */
    private String svgHBar(List<String> labels, List<Double> values, String title, String xLabel, String color, int svgW) {
        if (labels.isEmpty()) return "";
        int n = labels.size();
        int labelW = 200;
        int hBarH = 22, hGap = 6;
        int svgH = MT + n * (hBarH + hGap) + MB + 10;
        int pw = svgW - labelW - MR;
        int ph = n * (hBarH + hGap);

        double maxVal = values.stream().mapToDouble(Math::abs).max().orElse(1);
        double[] xs = niceScale(0, maxVal, 5);
        double xMax = xs[1], xStep = xs[2];
        if (xMax == 0) xMax = 1;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" style=\"width:100%%;max-width:%dpx\">\n", svgW, svgH, svgW));

        if (!title.isEmpty()) sb.append(svgText(svgW / 2.0, 20, title, "middle", 13, true));

        // Grid
        for (double tick = xStep; tick <= xMax; tick += xStep) {
            double x = labelW + tick / xMax * pw;
            sb.append(String.format("<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#e8e8e8\" stroke-width=\"1\"/>\n", x, MT, x, MT + ph));
        }

        // Bars
        for (int i = 0; i < n; i++) {
            double v = values.get(i);
            double barPx = Math.abs(v) / xMax * pw;
            double y = MT + i * (hBarH + hGap);
            sb.append(String.format("<rect x=\"%d\" y=\"%.1f\" width=\"%.1f\" height=\"%d\" fill=\"%s\" rx=\"2\"/>\n",
                    labelW, y, barPx, hBarH, color));
            sb.append(svgText(labelW - 5, y + hBarH / 2.0 + 4, labels.get(i), "end", 10, false));
            sb.append(svgText(labelW + barPx + 4, y + hBarH / 2.0 + 4, String.format("%.2f", v), "start", 9, false));
        }

        // X axis
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", labelW, MT + ph, labelW + pw, MT + ph));

        // X ticks
        for (double tick = 0; tick <= xMax + xStep * 0.01; tick += xStep) {
            double x = labelW + tick / xMax * pw;
            sb.append(svgText(x, MT + ph + 15, fmtTick(tick), "middle", 9, false));
        }

        sb.append(svgText(labelW + pw / 2.0, svgH - 5, xLabel, "middle", 11, false));
        sb.append("</svg>\n");
        return sb.toString();
    }

    /** Horizontal bar chart with per-bar colors (for percolator features with +/- values) */
    private String svgHBarColored(List<String> labels, List<Double> values, List<String> colors, String title, String xLabel, int svgW) {
        if (labels.isEmpty()) return "";
        int n = labels.size();
        int labelW = 180;
        int hBarH = 20, hGap = 5;
        int svgH = MT + n * (hBarH + hGap) + MB + 10;
        int pw = svgW - labelW - MR;
        int ph = n * (hBarH + hGap);

        double absMax = values.stream().mapToDouble(Math::abs).max().orElse(1);
        double[] xs = niceScale(-absMax, absMax, 5);
        double xLo = xs[0], xHi = xs[1], xStep = xs[2];
        double xRange = xHi - xLo;
        if (xRange == 0) xRange = 1;

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" style=\"width:100%%;max-width:%dpx\">\n", svgW, svgH, svgW));
        if (!title.isEmpty()) sb.append(svgText(svgW / 2.0, 20, title, "middle", 13, true));

        // Zero line
        double zeroX = labelW + (0 - xLo) / xRange * pw;
        sb.append(String.format("<line x1=\"%.1f\" y1=\"%d\" x2=\"%.1f\" y2=\"%d\" stroke=\"#999\" stroke-width=\"1\" stroke-dasharray=\"4,2\"/>\n", zeroX, MT, zeroX, MT + ph));

        // Bars
        for (int i = 0; i < n; i++) {
            double v = values.get(i);
            double vx = labelW + (v - xLo) / xRange * pw;
            double y = MT + i * (hBarH + hGap);
            double bx = Math.min(zeroX, vx);
            double bw = Math.abs(vx - zeroX);
            sb.append(String.format("<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%d\" fill=\"%s\" rx=\"2\"/>\n",
                    bx, y, bw, hBarH, colors.get(i)));
            sb.append(svgText(labelW - 5, y + hBarH / 2.0 + 4, labels.get(i), "end", 9, false));
            sb.append(svgText(vx + (v >= 0 ? 4 : -4), y + hBarH / 2.0 + 4, String.format("%.4f", v), v >= 0 ? "start" : "end", 8, false));
        }

        // X axis
        sb.append(String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#333\" stroke-width=\"1\"/>\n", labelW, MT + ph, labelW + pw, MT + ph));

        // X ticks
        for (double tick = xLo; tick <= xHi + xStep * 0.01; tick += xStep) {
            double x = labelW + (tick - xLo) / xRange * pw;
            sb.append(svgText(x, MT + ph + 15, fmtTick(tick), "middle", 9, false));
        }

        sb.append(svgText(labelW + pw / 2.0, svgH - 5, xLabel, "middle", 11, false));
        sb.append("</svg>\n");
        return sb.toString();
    }

    // ============================================================
    // SVG helpers
    // ============================================================

    private static double[] niceScale(double dataMin, double dataMax, int maxTicks) {
        if (dataMin == dataMax) { dataMax = dataMin + 1; }
        double range = dataMax - dataMin;
        double roughStep = range / maxTicks;
        double mag = Math.pow(10, Math.floor(Math.log10(roughStep)));
        double residual = roughStep / mag;
        double niceStep;
        if (residual <= 1.5) niceStep = mag;
        else if (residual <= 3) niceStep = 2 * mag;
        else if (residual <= 7) niceStep = 5 * mag;
        else niceStep = 10 * mag;
        double niceMin = Math.floor(dataMin / niceStep) * niceStep;
        double niceMax = Math.ceil(dataMax / niceStep) * niceStep;
        return new double[]{niceMin, niceMax, niceStep};
    }

    private static String fmtTick(double v) {
        if (v == Math.floor(v) && !Double.isInfinite(v) && Math.abs(v) < 1e15) return String.valueOf((long) v);
        if (Math.abs(v) >= 100) return String.format("%.0f", v);
        if (Math.abs(v) >= 1) return String.format("%.1f", v);
        return String.format("%.2f", v);
    }

    private static String svgText(double x, double y, String text, String anchor, int fontSize, boolean bold) {
        return String.format("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"%s\" font-size=\"%d\"%s>%s</text>\n",
                x, y, anchor, fontSize, bold ? " font-weight=\"bold\"" : "", esc(text));
    }

    private static String svgRotText(double x, double y, String text, int fontSize) {
        return String.format("<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" font-size=\"%d\" transform=\"rotate(-90,%.1f,%.1f)\">%s</text>\n",
                x, y, fontSize, x, y, esc(text));
    }

    private static String esc(String s) {
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\"", "&quot;");
    }

    // ============================================================
    // Utility methods
    // ============================================================

    private String imageToBase64(String path) {
        try {
            byte[] bytes = Files.readAllBytes(Path.of(path));
            return Base64.getEncoder().encodeToString(bytes);
        } catch (Exception e) { return null; }
    }

    private boolean isNoExpNoBio() { return hasAllMissing("Experiment") && hasAllMissing("Bioreplicate"); }
    private boolean isBioNull() { return hasAllMissing("Bioreplicate") && !hasAllMissing("Experiment"); }
    private boolean isExpNull() { return hasAllMissing("Experiment") && !hasAllMissing("Bioreplicate"); }

    private boolean hasAllMissing(String colName) {
        for (int i = 0; i < manifestData.rowCount(); i++) {
            String val = getString(manifestData, colName, i);
            if (val != null && !val.trim().isEmpty()) return false;
        }
        return true;
    }

    private List<String> uniqueExperiments() {
        Set<String> seen = new LinkedHashSet<>();
        for (int i = 0; i < manifestData.rowCount(); i++) {
            String v = getString(manifestData, "Experiment", i);
            if (v != null && !v.trim().isEmpty()) seen.add(v);
        }
        return new ArrayList<>(seen);
    }

    private List<String> uniqueBioreplicates() {
        Set<String> seen = new LinkedHashSet<>();
        for (int i = 0; i < manifestData.rowCount(); i++) {
            String v = getString(manifestData, "Bioreplicate", i);
            if (v != null && !v.trim().isEmpty()) seen.add(v);
        }
        return new ArrayList<>(seen);
    }

    private int countLines(String filename, String exp) {
        try {
            Table t = Table.read().csv(CsvReadOptions.builder(filename).separator('\t').build());
            convertTextToString(t);
            if (!exp.isEmpty() && t.columnNames().contains("Peptide Length")
                    && t.columnNames().contains("Charge") && t.columnNames().contains("Number of Missed Cleavages")) {
                // Deduplicate by (Spectrum, Modified Peptide, Peptide Length, Charge, Missed Cleavages)
                // to match the Python version's drop_duplicates behavior
                boolean canDedup = t.columnNames().contains("Spectrum") && t.columnNames().contains("Modified Peptide");
                Set<String> seen = canDedup ? new HashSet<>() : null;
                for (int i = 0; i < t.rowCount(); i++) {
                    if (canDedup) {
                        String key = t.column("Spectrum").getString(i) + "\t"
                                + t.column("Modified Peptide").getString(i) + "\t"
                                + t.column("Peptide Length").getString(i) + "\t"
                                + t.column("Charge").getString(i) + "\t"
                                + t.column("Number of Missed Cleavages").getString(i);
                        if (!seen.add(key)) continue;
                    }
                    distributionData.doubleColumn("Peptide Length").append(toDouble(t, "Peptide Length", i));
                    distributionData.doubleColumn("Charge").append(toDouble(t, "Charge", i));
                    distributionData.doubleColumn("Number of Missed Cleavages").append(toDouble(t, "Number of Missed Cleavages", i));
                    distributionData.stringColumn("Exp").append(exp);
                }
            }
            return t.rowCount();
        } catch (Exception e) { return 0; }
    }

    private void addIdRow(String experiment, double psm, double peptides, double proteins) {
        idNums.stringColumn("Experiment").append(experiment);
        idNums.doubleColumn("PSM").append(psm);
        idNums.doubleColumn("Peptides").append(peptides);
        idNums.doubleColumn("Proteins").append(proteins);
    }

    private void appendTable(Table target, Table source) {
        for (String colName : target.columnNames()) {
            Column<?> tgtCol = target.column(colName);
            if (source.columnNames().contains(colName)) {
                Column<?> srcCol = source.column(colName);
                for (int i = 0; i < srcCol.size(); i++) {
                    if (tgtCol instanceof StringColumn) ((StringColumn) tgtCol).append(srcCol.getString(i));
                    else if (tgtCol instanceof DoubleColumn) {
                        if (srcCol instanceof DoubleColumn) ((DoubleColumn) tgtCol).append(((DoubleColumn) srcCol).get(i));
                        else ((DoubleColumn) tgtCol).appendMissing();
                    }
                }
            } else {
                for (int i = 0; i < source.rowCount(); i++) {
                    if (tgtCol instanceof DoubleColumn) ((DoubleColumn) tgtCol).appendMissing();
                    else if (tgtCol instanceof StringColumn) ((StringColumn) tgtCol).append("");
                }
            }
        }
    }

    private Table tableFromMaps(String name, List<Map<String, Object>> rows) {
        StringColumn rc = StringColumn.create("Run"), cc = StringColumn.create("Calib");
        DoubleColumn vc = DoubleColumn.create("Value");
        for (Map<String, Object> r : rows) { rc.append((String) r.get("Run")); cc.append((String) r.get("Calib")); vc.append((Double) r.get("Value")); }
        return Table.create(name, rc, cc, vc);
    }

    private double toDouble(Table t, String colName, int row) {
        try {
            Column<?> col = t.column(colName);
            if (col instanceof DoubleColumn) return ((DoubleColumn) col).get(row);
            if (col instanceof IntColumn) return ((IntColumn) col).get(row);
            return Double.parseDouble(col.getString(row));
        } catch (Exception e) { return Double.NaN; }
    }

    private double[] toDoubleArrayNoNaN(DoubleColumn col) {
        List<Double> v = new ArrayList<>();
        for (int i = 0; i < col.size(); i++) if (!col.isMissing(i) && !Double.isNaN(col.get(i))) v.add(col.get(i));
        return v.stream().mapToDouble(Double::doubleValue).toArray();
    }

    /** Group idNums by Experiment, taking the first non-NaN value per metric column. */
    private Table groupIdNums() {
        LinkedHashMap<String, double[]> map = new LinkedHashMap<>();
        for (int i = 0; i < idNums.rowCount(); i++) {
            String exp = idNums.stringColumn("Experiment").get(i);
            double[] vals = map.computeIfAbsent(exp, k -> new double[]{Double.NaN, Double.NaN, Double.NaN});
            for (int c = 0; c < 3; c++) {
                String col = new String[]{"PSM", "Peptides", "Proteins"}[c];
                if (!idNums.columnNames().contains(col)) continue;
                DoubleColumn dc = idNums.doubleColumn(col);
                if (!dc.isMissing(i) && !Double.isNaN(dc.get(i)) && Double.isNaN(vals[c]))
                    vals[c] = dc.get(i);
            }
        }
        StringColumn expCol = StringColumn.create("Experiment");
        DoubleColumn psmCol = DoubleColumn.create("PSM");
        DoubleColumn pepCol = DoubleColumn.create("Peptides");
        DoubleColumn protCol = DoubleColumn.create("Proteins");
        for (Map.Entry<String, double[]> e : map.entrySet()) {
            expCol.append(e.getKey());
            double[] v = e.getValue();
            if (Double.isNaN(v[0])) psmCol.appendMissing(); else psmCol.append(v[0]);
            if (Double.isNaN(v[1])) pepCol.appendMissing(); else pepCol.append(v[1]);
            if (Double.isNaN(v[2])) protCol.appendMissing(); else protCol.append(v[2]);
        }
        return Table.create("grouped", expCol, psmCol, pepCol, protCol);
    }

    private boolean hasCol(Table t, String name) { return t.columnNames().contains(name); }

    // ============================================================
    // Main
    // ============================================================

    public static void main(String[] args) {
        String resultsPath = null;
        for (int i = 0; i < args.length; i++)
            if (("-r".equals(args[i]) || "--results_path".equals(args[i])) && i + 1 < args.length) { resultsPath = args[++i]; }

        if (resultsPath == null) { System.err.println("Usage: java -jar fragsummarizer.jar -r <results_path>"); System.exit(1); }
        if (!Files.isDirectory(Path.of(resultsPath))) { System.err.println("Error: Not a directory: " + resultsPath); System.exit(1); }

        new FragSummarizer(resultsPath).generateReport();
    }
}
