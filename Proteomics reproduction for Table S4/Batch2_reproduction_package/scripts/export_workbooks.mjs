import { readFile } from "node:fs/promises";
import { createRequire } from "node:module";
import path from "node:path";
import { pathToFileURL } from "node:url";

const [directory] = process.argv.slice(2);
if (!directory) throw new Error("Supply a directory containing the generated CSV files.");
const modules = process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES;
const specifier = modules
  ? pathToFileURL(createRequire(path.join(modules, "..", "artifact-loader.cjs")).resolve("@oai/artifact-tool")).href
  : "@oai/artifact-tool";
const { SpreadsheetFile, Workbook } = await import(specifier);

for (const [stem, name] of [
  ["Batch2_PPM_full", "Batch2 PPM"],
  ["TableS4_Protein_reproduced", "dnaA&hns_Proteins"],
]) {
  const csv = (await readFile(path.join(directory, stem + ".csv"), "utf8")).replace(/^\uFEFF/, "");
  const workbook = await Workbook.fromCSV(csv, { sheetName: name });
  const sheet = workbook.worksheets.getItem(name);
  const range = sheet.getUsedRange();
  const values = range.values;
  const numeric = values[0].map((label) => /^(dnaA|hns)upshift_Rep|^Is_(division|dcw)_genes$/.test(label));
  for (const row of values.slice(1)) {
    numeric.forEach((isNumeric, column) => {
      if (row[column] === "" || row[column] == null) row[column] = null;
      else if (isNumeric && row[column] !== "NA") row[column] = Number(row[column]);
    });
  }
  range.values = values;
  sheet.freezePanes.freezeRows(1);
  const header = sheet.getRangeByIndexes(0, 0, 1, values[0].length);
  header.format.fill = "#163958";
  header.format.font.color = "#FFFFFF";
  header.format.font.bold = true;
  header.format.wrapText = true;
  header.format.rowHeight = 40;
  sheet.getRange("A:A").format.columnWidth = 23;
  sheet.getRange("B:B").format.columnWidth = 18;
  sheet.getRange("C:F").format.columnWidth = 15;
  sheet.getRange("G:G").format.columnWidth = 39;
  sheet.getRange("H:BY").format.columnWidth = 17;
  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(path.join(directory, stem + ".xlsx"));
}
