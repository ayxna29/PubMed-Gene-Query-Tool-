# PubMed Gene Query Tool

A Shiny application for batch querying PubMed across a gene list and multiple keyword sets, returning the top N studies per gene with clickable links.

---

## What It Does

You upload a list of genes, define one or more keyword search sets (e.g. `ovarian AND cancer`, `drug resistance AND therapy`), and the app queries PubMed for every gene–keyword combination. Results are displayed in an interactive table where each cell shows a study title linked directly to its PubMed page. You can export to CSV or Excel.

---

## Requirements

### R Packages
```r
install.packages(c("shiny", "rentrez", "DT", "openxlsx"))
```

- **shiny** — web application framework
- **rentrez** — NCBI Entrez API interface
- **DT** — interactive DataTables for Shiny
- **openxlsx** — Excel file export (optional; only needed for `.xlsx` download)

---

## How to Run

```r
library(shiny)
runApp("pubmed_gene_query.R")
```

Or open the file in RStudio and click **Run App**.

---

## Input File Format

Upload a `.txt` or `.csv` file with gene names in the **first column**, one per row.

**CSV example:**
```
Gene
BRCA1
TP53
EGFR
PTEN
```

**TXT example (tab-separated or one per line):**
```
Gene
BRCA1
TP53
EGFR
```

A header row is expected. If your file has no header, the first gene will be skipped.

---

## Usage

1. **Upload** your gene list using the file input
2. **Configure** search options:
   - Sort by relevance or most recent
   - Set how many top results to retrieve per gene (1–10)
   - Optionally restrict to a publication year range
3. **Add keyword sets** — each set becomes one PubMed search combined with every gene as `[Gene] AND ([your keywords])`. Click **+ Add Keyword Set** to add more. Use the ✕ button to remove any set.
4. **Set the delay** between queries (default 0.35 s — keep at or above this to respect NCBI rate limits)
5. **Click Run PubMed Query** — a progress bar tracks each gene across each keyword set
6. **View results** in the Results tab; use the Summary tab for a high-level overview
7. **Download** as CSV or Excel

---

## Output

### In-App Table
Each `_Top1`, `_Top2`, `_Top3` column (per keyword set) shows the study title as a clickable hyperlink to its PubMed page. Hovering over a link shows year, journal, and first author in a tooltip.

### Column Naming
Columns are named using your actual keyword text, sanitised for use as column names. For example, a keyword set `breast AND cancer` produces columns:

```
breast_AND_cancer_Top1
breast_AND_cancer_Top2
breast_AND_cancer_Top3
breast_AND_cancer_TotalHits
```

### CSV Download
Each `_TopN` column contains the plain study title (no URL — CSV cannot embed clickable links).

### Excel Download
Each `_TopN` column contains a real `=HYPERLINK()` formula — the study title appears as blue underlined text and clicking it opens the PubMed page directly.

---

## Tabs

| Tab | Contents |
|---|---|
| Gene Preview | Shows how many genes were loaded and a preview of the first 30 |
| Results | Interactive filterable table with all query results |
| Summary | Bar chart and table showing how many genes had ≥1 hit per keyword set |

---

## NCBI Rate Limits

NCBI allows up to 3 requests/second without an API key, and 10/second with one. The default delay of 0.35 s keeps usage safely within the unauthenticated limit. For large gene lists (500+), consider registering for a free NCBI API key and setting the delay lower:

```r
rentrez::set_entrez_key("YOUR_API_KEY")
```

See [NCBI API Key documentation](https://www.ncbi.nlm.nih.gov/account/) for details.

---

## Tips

- Use PubMed Boolean syntax in keyword sets: `AND`, `OR`, `NOT`, field tags like `[MeSH Terms]`, `[tiab]`
- Example keyword sets:
  - `cancer AND therapy`
  - `"drug resistance" AND ("ovarian cancer"[MeSH Terms])`
  - `expression AND "RNA-seq"[tiab]`
- The **Total Hits** column shows the total number of PubMed results for that gene–keyword combination, not just the top N retrieved
- The debug line beneath the keyword boxes confirms what text is registered before you click Run

---

## Acknowledgements
**Tianyang Li** — authored the original R script that forms the foundation of 
  this application, including the core PubMed query function and initial Shiny app skeleton with file upload, 
  keyword inputs, progress tracking, and CSV export.
Extended by Ayana Singh 
