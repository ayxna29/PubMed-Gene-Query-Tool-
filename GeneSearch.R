library(shiny)
library(rentrez)
library(DT)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Gene symbol → full name (local Bioconductor DB, no API needed) ────────────
# Requires: BiocManager::install("org.Hs.eg.db")
# Returns named character vector: names = original symbols, values = full names.
# Falls back to the original symbol for any gene not found.
convert_gene_symbols <- function(symbols) {
  result <- setNames(symbols, symbols)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("org.Hs.eg.db not installed — using raw symbols.")
    return(result)
  }
  tryCatch({
    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = symbols,
      columns = "GENENAME",
      keytype = "SYMBOL"
    )
    map <- map[!duplicated(map$SYMBOL), ]
    found <- !is.na(map$GENENAME) & nchar(trimws(map$GENENAME)) > 0
    result[map$SYMBOL[found]] <- map$GENENAME[found]
  }, error = function(e) message("Gene name conversion error: ", e$message))
  result
}

# ── helpers ───────────────────────────────────────────────────────────────────

make_link <- function(pmid, title, year = NA, journal = NA, authors = NA) {
  if (is.na(pmid) || pmid == "" || is.na(title) || title == "") return(NA_character_)
  url   <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/")
  title <- gsub("'", "&#39;", title)
  tip_parts <- c(
    if (!is.na(year)    && nchar(year)    > 0) paste0("Year: ",    year),
    if (!is.na(journal) && nchar(journal) > 0) paste0("Journal: ", journal),
    if (!is.na(authors) && nchar(authors) > 0) paste0("Author: ",  authors)
  )
  paste0("<a href='", url, "' target='_blank' title='", paste(tip_parts, collapse=" | "), "'>", title, "</a>")
}

get_top_pubmed_results <- function(query, n = 3, sort_by = "relevance") {
  empty <- list(ids=character(0), titles=character(0), links=character(0),
                total_count=0L, years=character(0), journals=character(0), urls=character(0))
  sort_arg <- if (sort_by == "most_recent") "pub_date" else "relevance"
  res <- tryCatch(entrez_search(db="pubmed", term=query, retmax=n, sort=sort_arg),
                  error = function(e) { message("entrez_search error: ", e$message); NULL })
  if (is.null(res) || length(res$ids) == 0) return(empty)
  total_count <- res$count
  ids <- res$ids[seq_len(min(n, length(res$ids)))]
  smm <- tryCatch(entrez_summary(db="pubmed", id=ids),
                  error = function(e) { message("entrez_summary error: ", e$message); NULL })
  extract_field <- function(item, field) {
    val <- tryCatch(item[[field]], error = function(e) NA_character_)
    if (is.null(val) || length(val) == 0) NA_character_ else as.character(val[1])
  }
  get_authors <- function(item) {
    au <- tryCatch(item$authors, error = function(e) NULL)
    if (!is.null(au) && is.data.frame(au) && nrow(au) > 0) paste0(au$name[1], " et al.")
    else NA_character_
  }
  if (is.null(smm)) {
    titles <- years <- journals <- authors <- rep(NA_character_, length(ids))
  } else if (length(ids) == 1) {
    titles   <- extract_field(smm, "title")
    years    <- extract_field(smm, "pubdate")
    journals <- extract_field(smm, "fulljournalname")
    authors  <- get_authors(smm)
  } else {
    titles   <- vapply(ids, function(id) { item <- tryCatch(smm[[id]], error=function(e) NULL); if (is.null(item)) NA_character_ else extract_field(item, "title")          }, character(1), USE.NAMES=FALSE)
    years    <- vapply(ids, function(id) { item <- tryCatch(smm[[id]], error=function(e) NULL); if (is.null(item)) NA_character_ else extract_field(item, "pubdate")         }, character(1), USE.NAMES=FALSE)
    journals <- vapply(ids, function(id) { item <- tryCatch(smm[[id]], error=function(e) NULL); if (is.null(item)) NA_character_ else extract_field(item, "fulljournalname")  }, character(1), USE.NAMES=FALSE)
    authors  <- vapply(ids, function(id) { item <- tryCatch(smm[[id]], error=function(e) NULL); if (is.null(item)) NA_character_ else get_authors(item)                     }, character(1), USE.NAMES=FALSE)
  }
  years <- sub("^(\\d{4}).*", "\\1", years)
  urls  <- paste0("https://pubmed.ncbi.nlm.nih.gov/", ids, "/")
  links <- mapply(make_link, ids, titles, years, journals, authors, USE.NAMES=FALSE)
  list(ids=ids, titles=titles, links=links, urls=urls,
       total_count=total_count, years=years, journals=journals)
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background:#f7f8fa; font-family:'Segoe UI',sans-serif; }
    .main-title { color:#2c3e50; margin-bottom:4px; }
    .subtitle   { color:#7f8c8d; margin-top:0; margin-bottom:20px; }
    .well       { background:#fff; border:1px solid #dee2e6; }
    .btn-primary { background:#3498db; border:none; }
    .btn-primary:hover { background:#2980b9; }
    #status_box { margin-top:10px; padding:8px 12px; border-radius:6px;
                  background:#eaf4fb; border:1px solid #aed6f1;
                  font-size:13px; color:#1a5276; }
    .kw-block   { background:#f0f4f8; border-radius:6px;
                  padding:8px 10px; margin-bottom:6px; position:relative; }
    .kw-label   { font-size:12px; font-weight:600; color:#2c3e50; margin-bottom:4px; }
    .remove-btn { position:absolute; top:8px; right:8px; padding:1px 7px;
                  font-size:12px; background:#e74c3c; color:#fff; border:none;
                  border-radius:4px; cursor:pointer; line-height:1.6; }
    a { color:#2980b9; }
    .debug-box  { font-size:11px; color:#888; font-family:monospace; margin-top:4px; }
  "))),
  
  titlePanel(div(
    h2("PubMed Gene Query Tool", class="main-title"),
    p("Upload a gene list, configure keyword searches, and retrieve top PubMed studies.", class="subtitle")
  )),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      fileInput("gene_file", "Upload Gene List (.txt or .csv)",
                accept=c("text/plain","text/csv",".txt",".csv")),
      helpText("One gene per row; gene names in the first column."),
      
      checkboxInput("convert_symbols",
                    "Also search gene full names (e.g. ACKR3 + 'atypical chemokine receptor 3')",
                    value = TRUE),
      helpText("Searches as: (SYMBOL OR \"full name\") AND keywords."),
      
      hr(),
      h5("Search Options"),
      selectInput("sort_by", "Sort results by",
                  choices=c("Relevance"="relevance","Most Recent"="most_recent"), selected="relevance"),
      numericInput("top_n", "Top N results per gene per keyword", value=3, min=1, max=10, step=1),
      checkboxInput("include_total", "Add 'Total Hits' column", value=TRUE),
      
      hr(),
      h5("Date Filter (optional)"),
      fluidRow(
        column(6, numericInput("year_from", "From year", value=NA, min=1900, max=2100, step=1)),
        column(6, numericInput("year_to",   "To year",   value=NA, min=1900, max=2100, step=1))
      ),
      helpText("Leave blank for no date restriction."),
      
      hr(),
      actionButton("add_combo", "+ Add keyword set", class="btn-primary btn-block"),
      helpText("Each box: [Gene] AND [your keywords]."),
      br(),
      uiOutput("keyword_ui"),
      div(class="debug-box", textOutput("active_kw_debug")),
      
      hr(),
      numericInput("delay", "Delay between queries (seconds)", value=0.35, min=0.1, max=5, step=0.05),
      helpText("Keep >= 0.34 s to respect NCBI rate limits."),
      
      hr(),
      actionButton("run_btn", "Run PubMed Query", icon=icon("search"), class="btn-primary btn-block"),
      br(),
      div(id="status_box", textOutput("status_text")),
      
      hr(),
      downloadButton("download_btn",      "Download CSV",            class="btn-block"),
      br(),
      downloadButton("download_xlsx_btn", "Download Excel (.xlsx)",  class="btn-block")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Gene Preview", br(), verbatimTextOutput("gene_preview")),
        tabPanel("Results",
                 br(),
                 uiOutput("col_visibility_ui"),
                 br(),
                 DTOutput("results_table")),
        tabPanel("Summary",
                 br(),
                 plotOutput("hit_plot", height="400px"),
                 br(),
                 DTOutput("summary_table"))
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  # ── Keyword slots ──────────────────────────────────────────────────────────
  kw_slots     <- reactiveVal(list(list(slot_id="kw_1", label="Keyword Set 1")))
  next_slot_id <- reactiveVal(2L)
  
  observeEvent(input$add_combo, {
    sid <- paste0("kw_", next_slot_id())
    kw_slots(c(kw_slots(), list(list(slot_id=sid, label=paste("Keyword Set", next_slot_id())))))
    next_slot_id(next_slot_id() + 1L)
  })
  observeEvent(input$remove_kw, {
    kw_slots(Filter(function(s) s$slot_id != input$remove_kw, kw_slots()))
  })
  
  output$keyword_ui <- renderUI({
    tagList(lapply(kw_slots(), function(s) {
      sid <- s$slot_id
      div(class="kw-block",
          tags$button("x", class="remove-btn",
                      onclick=sprintf("Shiny.setInputValue('remove_kw','%s',{priority:'event'})", sid)),
          div(class="kw-label", s$label),
          textInput(paste0("text_", sid), label=NULL, placeholder="e.g. ovarian AND cancer AND therapy")
      )
    }))
  })
  
  read_active_keywords <- function() {
    result <- list()
    for (s in kw_slots()) {
      val <- trimws(input[[paste0("text_", s$slot_id)]] %||% "")
      if (nchar(val) > 0) result[[s$slot_id]] <- list(text=val, label=s$label)
    }
    result
  }
  
  output$active_kw_debug <- renderText({
    vals <- vapply(kw_slots(), function(s) {
      v <- trimws(input[[paste0("text_", s$slot_id)]] %||% "")
      if (nchar(v) > 0) paste0(s$label, ": ", substr(v,1,30)) else paste0(s$label, ": (empty)")
    }, character(1))
    paste("Active:", paste(vals, collapse=" | "))
  })
  
  # ── Gene list ──────────────────────────────────────────────────────────────
  gene_list <- reactive({
    req(input$gene_file)
    ext <- tolower(tools::file_ext(input$gene_file$name))
    df  <- if (ext == "csv") {
      read.csv(input$gene_file$datapath, header=TRUE, stringsAsFactors=FALSE)
    } else {
      tryCatch(
        read.table(input$gene_file$datapath, header=TRUE, stringsAsFactors=FALSE, sep="\t"),
        error = function(e)
          data.frame(Gene=readLines(input$gene_file$datapath), stringsAsFactors=FALSE)
      )
    }
    genes <- trimws(as.character(df[[1]]))
    genes[nchar(genes) > 0]
  })
  
  # ── Symbol → full name map ─────────────────────────────────────────────────
  symbol_map <- reactive({
    genes <- gene_list()
    if (!isTRUE(input$convert_symbols)) return(setNames(rep(NA_character_, length(genes)), genes))
    convert_gene_symbols(genes)
  })
  
  output$gene_preview <- renderPrint({
    genes <- gene_list()
    sm    <- symbol_map()
    cat(sprintf("Total genes loaded: %d\n\n", length(genes)))
    preview <- genes[seq_len(min(30, length(genes)))]
    if (isTRUE(input$convert_symbols)) {
      for (g in preview) {
        nm <- sm[[g]]
        if (!is.na(nm) && g != nm) cat(sprintf("  %-15s + \"%s\"\n", g, nm))
        else                        cat(sprintf("  %s  (symbol only)\n", g))
      }
    } else {
      cat(paste(preview, collapse="\n"))
    }
    if (length(genes) > 30) cat(sprintf("\n... and %d more", length(genes) - 30))
  })
  
  # ── Date filter ────────────────────────────────────────────────────────────
  date_filter <- reactive({
    yf <- input$year_from; yt <- input$year_to
    if (!is.na(yf) && !is.na(yt)) sprintf(" AND %d:%d[pdat]", yf, yt)
    else if (!is.na(yf))          sprintf(" AND %d:3000[pdat]", yf)
    else if (!is.na(yt))          sprintf(" AND 1900:%d[pdat]", yt)
    else ""
  })
  
  run_keywords <- reactiveVal(NULL)
  
  make_prefix <- function(kw_text) {
    p <- gsub("[^A-Za-z0-9]+", "_", kw_text)
    p <- gsub("^_|_$", "", p)
    substr(p, 1, 40)
  }
  
  # ── Main query ─────────────────────────────────────────────────────────────
  results <- eventReactive(input$run_btn, {
    genes     <- gene_list()
    sm        <- isolate(symbol_map())   # named vec: symbol -> full name (or NA)
    kws       <- read_active_keywords()
    topn      <- isolate(input$top_n)
    delay     <- isolate(input$delay)
    df_fil    <- isolate(date_filter())
    srt       <- isolate(input$sort_by)
    inc_total <- isolate(input$include_total)
    use_conv  <- isolate(input$convert_symbols)
    
    validate(
      need(length(genes) > 0, "Please upload a gene list first."),
      need(length(kws)   > 0, "Please enter at least one keyword set.")
    )
    run_keywords(kws)
    
    gene_names <- unname(sm[genes])   # full names, NA if not converted
    out <- data.frame(Gene=genes, Gene_Name=ifelse(is.na(gene_names), genes, gene_names),
                      stringsAsFactors=FALSE)
    
    for (j in seq_along(names(kws))) {
      sid        <- names(kws)[j]
      kw_text    <- kws[[sid]]$text
      col_prefix <- make_prefix(kw_text)
      
      # Column order: Top1_Title, Top1_URL, Top2_Title, Top2_URL, ...
      slot_cols <- unlist(lapply(seq_len(topn), function(k)
        c(paste0(col_prefix, "_Top", k, "_Title"),
          paste0(col_prefix, "_Top", k, "_URL"),
          paste0(col_prefix, "_Top", k, "_Link"))   # Link = HTML, for in-app table only
      ))
      hit_col <- paste0(col_prefix, "_TotalHits")
      
      out[slot_cols] <- NA_character_
      if (inc_total) out[[hit_col]] <- NA_integer_
      
      withProgress(
        message = sprintf("Keyword set %d/%d: %s", j, length(kws), substr(kw_text, 1, 30)),
        value = 0, {
          for (i in seq_along(genes)) {
            sym      <- genes[i]
            fullname <- gene_names[i]
            
            # Build query: (SYMBOL OR "full name") AND (keywords) — or just SYMBOL if no conversion
            gene_part <- if (use_conv && !is.na(fullname) && fullname != sym) {
              sprintf('(%s OR "%s")', sym, fullname)
            } else {
              sym
            }
            q   <- paste0(gene_part, " AND (", kw_text, ")", df_fil)
            top <- get_top_pubmed_results(q, n=topn, sort_by=srt)
            
            for (k in seq_len(topn)) {
              tc <- paste0(col_prefix, "_Top", k, "_Title")
              uc <- paste0(col_prefix, "_Top", k, "_URL")
              lc <- paste0(col_prefix, "_Top", k, "_Link")
              if (length(top$titles) >= k && !is.na(top$titles[k])) {
                out[i, tc] <- top$titles[k]
                out[i, uc] <- if (length(top$urls)  >= k) top$urls[k]  else NA_character_
                out[i, lc] <- if (length(top$links) >= k) top$links[k] else NA_character_
              }
            }
            if (inc_total) out[i, hit_col] <- top$total_count
            
            incProgress(1/length(genes), detail=sprintf("%d/%d — %s", i, length(genes), sym))
            Sys.sleep(delay)
          }
        }
      )
    }
    out
  })
  
  output$status_text <- renderText({
    df <- results(); kws <- run_keywords()
    n_populated <- sum(!is.na(df[, grep("_Title$", names(df)), drop=FALSE]))
    sprintf("Done — %d genes x %d keyword set(s) | %d results found.",
            nrow(df), length(kws), n_populated)
  })
  
  # ── Column visibility ──────────────────────────────────────────────────────
  output$col_visibility_ui <- renderUI({
    kws <- run_keywords()
    if (is.null(kws) || length(kws) <= 1) return(NULL)
    pfxs   <- vapply(kws, function(k) make_prefix(k$text), character(1))
    labels <- paste0("Set ", seq_along(kws), ": ",
                     substr(vapply(kws, `[[`, character(1), "text"), 1, 35))
    checkboxGroupInput("visible_searches", "Show keyword sets:",
                       choices=setNames(pfxs, labels), selected=pfxs, inline=TRUE)
  })
  
  filtered_results <- reactive({
    df  <- results(); vis <- input$visible_searches
    link_pat <- if (is.null(vis)) "_Link$" else
      paste(paste0("^", vis, "_.*_Link$"), collapse="|")
    hits_pat <- if (is.null(vis)) "_TotalHits$" else
      paste(paste0("^", vis, "_TotalHits$"), collapse="|")
    keep <- c("Gene", "Gene_Name",
              grep(link_pat, names(df), value=TRUE),
              grep(hits_pat, names(df), value=TRUE))
    df[, intersect(keep, names(df)), drop=FALSE]
  })
  
  output$results_table <- renderDT({
    df <- filtered_results()
    names(df) <- gsub("_Link$", "", names(df))
    datatable(df, rownames=FALSE, escape=FALSE, filter="top",
              options=list(pageLength=25, scrollX=TRUE,
                           columnDefs=list(list(width="300px",
                                                targets=grep("_Top\\d+$", names(df)) - 1L))))
  })
  
  # ── Summary tab ────────────────────────────────────────────────────────────
  summary_df <- reactive({
    df <- results(); kws <- run_keywords(); req(!is.null(kws))
    data.frame(
      Keyword_Set        = paste0("Set ", seq_along(kws)),
      Keywords           = vapply(kws, `[[`, character(1), "text"),
      Genes_With_Any_Hit = vapply(seq_along(kws), function(j) {
        cols <- grep(paste0("^", make_prefix(kws[[j]]$text), "_Top\\d+_Title$"), names(df), value=TRUE)
        if (!length(cols)) 0L else sum(rowSums(!is.na(df[cols, drop=FALSE])) > 0)
      }, integer(1)),
      Total_Links        = vapply(seq_along(kws), function(j) {
        cols <- grep(paste0("^", make_prefix(kws[[j]]$text), "_Top\\d+_Title$"), names(df), value=TRUE)
        if (!length(cols)) 0L else sum(!is.na(df[cols, drop=FALSE]))
      }, integer(1)),
      stringsAsFactors = FALSE
    )
  })
  
  output$summary_table <- renderDT({
    datatable(summary_df(), rownames=FALSE, options=list(dom="t"))
  })
  output$hit_plot <- renderPlot({
    sm <- summary_df()
    par(mar=c(5, 12, 3, 2))
    barplot(sm$Genes_With_Any_Hit,
            names.arg=paste0(sm$Keyword_Set, "\n(", substr(sm$Keywords, 1, 30), ")"),
            horiz=TRUE, las=1, col="#3498db", border=NA,
            main="Genes with >=1 PubMed hit per keyword set", xlab="Number of genes")
  })
  
  # ── CSV download: Title then URL interleaved, no HTML ─────────────────────
  output$download_btn <- downloadHandler(
    filename = function() paste0("pubmed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      df <- results()
      # Keep Gene, Gene_Name, TotalHits; interleave TopN_Title + TopN_URL per slot
      base_cols  <- c("Gene", "Gene_Name")
      title_cols <- grep("_Top\\d+_Title$", names(df), value=TRUE)
      url_cols   <- grep("_Top\\d+_URL$",   names(df), value=TRUE)
      hits_cols  <- grep("_TotalHits$",     names(df), value=TRUE)
      
      # Interleave: Top1_Title, Top1_URL, Top2_Title, Top2_URL, ...
      # They're already stored in interleaved order in 'df', just drop _Link cols
      ordered_link_cols <- unlist(lapply(seq_along(title_cols), function(k)
        c(title_cols[k], url_cols[k])
      ))
      
      keep <- c(base_cols, ordered_link_cols, hits_cols)
      write.csv(df[, intersect(keep, names(df)), drop=FALSE],
                file, row.names=FALSE, na="")
    }
  )
  
  # ── Excel download: title text visible + clickable hyperlink per cell ──────
  output$download_xlsx_btn <- downloadHandler(
    filename = function() paste0("pubmed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content  = function(file) {
      if (!requireNamespace("openxlsx", quietly=TRUE)) {
        showNotification("Install openxlsx: install.packages('openxlsx')", type="error")
        return()
      }
      library(openxlsx)
      df <- results()
      
      title_cols <- grep("_Top\\d+_Title$", names(df), value=TRUE)
      url_cols   <- grep("_Top\\d+_URL$",   names(df), value=TRUE)
      hits_cols  <- grep("_TotalHits$",     names(df), value=TRUE)
      base_cols  <- c("Gene", "Gene_Name")
      
      # Build display names: strip _Title suffix
      display_names <- sub("_Title$", "", title_cols)
      
      # Excel df: base + interleaved hyperlink columns + hits
      excel_df <- df[, base_cols, drop=FALSE]
      for (k in seq_along(title_cols)) {
        # Placeholder text = title; will be replaced by HYPERLINK formula below
        excel_df[[display_names[k]]] <- df[[title_cols[k]]]
      }
      for (hc in hits_cols) excel_df[[hc]] <- df[[hc]]
      
      wb <- createWorkbook()
      addWorksheet(wb, "PubMed Results")
      
      # Use writeData (not writeDataTable) so formulas aren't blocked by table formatting
      writeData(wb, 1, excel_df, startRow=1, startCol=1, headerStyle=
                  createStyle(fontColour="#FFFFFF", fgFill="#2c3e50", textDecoration="bold",
                              border="Bottom", borderColour="#aaaaaa"))
      
      # Freeze top row, auto-width
      freezePane(wb, 1, firstRow=TRUE)
      setColWidths(wb, 1, cols=seq_len(ncol(excel_df)), widths="auto")
      
      # Alternating row fill for readability
      for (row in seq_len(nrow(excel_df))) {
        fill_col <- if (row %% 2 == 0) "#f0f4f8" else "#ffffff"
        addStyle(wb, 1,
                 style    = createStyle(fgFill=fill_col),
                 rows     = row + 1L,
                 cols     = seq_len(ncol(excel_df)),
                 stack    = TRUE)
      }
      
      link_style <- createStyle(fontColour="#2980b9", textDecoration="underline")
      
      for (k in seq_along(title_cols)) {
        col_idx <- which(names(excel_df) == display_names[k])
        titles  <- df[[title_cols[k]]]
        urls    <- df[[url_cols[k]]]
        for (row in seq_len(nrow(excel_df))) {
          url   <- urls[row];   if (is.na(url)   || nchar(url)   == 0) next
          title <- titles[row]; if (is.na(title) || nchar(title) == 0) next
          # Write formula directly — no prior text in cell to conflict
          writeFormula(wb, 1,
                       x        = sprintf('HYPERLINK("%s","%s")', url, gsub('"', '""', title)),
                       startRow = row + 1L,
                       startCol = col_idx)
          addStyle(wb, 1, style=link_style, rows=row+1L, cols=col_idx, stack=TRUE)
        }
      }
      
      saveWorkbook(wb, file, overwrite=TRUE)
    }
  )
}

shinyApp(ui=ui, server=server)
