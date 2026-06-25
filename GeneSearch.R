library(shiny)
library(rentrez)
library(DT)

`%||%` <- function(a, b) if (!is.null(a)) a else b

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
  tip <- paste(tip_parts, collapse = " | ")
  
  paste0("<a href='", url, "' target='_blank' title='", tip, "'>", title, "</a>")
}

get_top_pubmed_results <- function(query, n = 3, sort_by = "relevance") {
  empty <- list(ids=character(0), titles=character(0), links=character(0),
                total_count=0L, years=character(0), journals=character(0))
  
  sort_arg <- if (sort_by == "most_recent") "pub_date" else "relevance"
  
  res <- tryCatch(
    entrez_search(db="pubmed", term=query, retmax=n, sort=sort_arg),
    error = function(e) { message("entrez_search error: ", e$message); NULL }
  )
  
  if (is.null(res) || length(res$ids) == 0) return(empty)
  
  total_count <- res$count
  ids <- res$ids[seq_len(min(n, length(res$ids)))]
  
  smm <- tryCatch(
    entrez_summary(db="pubmed", id=ids),
    error = function(e) { message("entrez_summary error: ", e$message); NULL }
  )
  
  extract_field <- function(item, field) {
    val <- tryCatch(item[[field]], error = function(e) NA_character_)
    if (is.null(val) || length(val) == 0) NA_character_ else as.character(val[1])
  }
  
  get_authors <- function(item) {
    au <- tryCatch(item$authors, error = function(e) NULL)
    if (!is.null(au) && is.data.frame(au) && nrow(au) > 0)
      paste0(au$name[1], " et al.")
    else
      NA_character_
  }
  
  if (is.null(smm)) {
    titles <- years <- journals <- authors <- rep(NA_character_, length(ids))
  } else if (length(ids) == 1) {
    # Single result — smm IS the summary object directly
    titles   <- extract_field(smm, "title")
    years    <- extract_field(smm, "pubdate")
    journals <- extract_field(smm, "fulljournalname")
    authors  <- get_authors(smm)
  } else {
    # Multiple results — smm is a named list keyed by PMID
    titles   <- vapply(ids, function(id) {
      item <- tryCatch(smm[[id]], error = function(e) NULL)
      if (is.null(item)) NA_character_ else extract_field(item, "title")
    }, character(1), USE.NAMES=FALSE)
    
    years    <- vapply(ids, function(id) {
      item <- tryCatch(smm[[id]], error = function(e) NULL)
      if (is.null(item)) NA_character_ else extract_field(item, "pubdate")
    }, character(1), USE.NAMES=FALSE)
    
    journals <- vapply(ids, function(id) {
      item <- tryCatch(smm[[id]], error = function(e) NULL)
      if (is.null(item)) NA_character_ else extract_field(item, "fulljournalname")
    }, character(1), USE.NAMES=FALSE)
    
    authors  <- vapply(ids, function(id) {
      item <- tryCatch(smm[[id]], error = function(e) NULL)
      if (is.null(item)) NA_character_ else get_authors(item)
    }, character(1), USE.NAMES=FALSE)
  }
  
  years <- sub("^(\\d{4}).*", "\\1", years)
  
  links <- mapply(make_link, ids, titles, years, journals, authors,
                  USE.NAMES=FALSE)
  
  list(ids=ids, titles=titles, links=links,
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
    .kw-label   { font-size:12px; font-weight:600; color:#2c3e50;
                  margin-bottom:4px; }
    .remove-btn { position:absolute; top:8px; right:8px; padding:1px 7px;
                  font-size:12px; background:#e74c3c; color:#fff; border:none;
                  border-radius:4px; cursor:pointer; line-height:1.6; }
    a { color:#2980b9; }
    .debug-box  { font-size:11px; color:#888; font-family:monospace;
                  margin-top:4px; }
  "))),
  
  titlePanel(div(
    h2("PubMed Gene Query Tool", class="main-title"),
    p("Upload a gene list, configure keyword searches, and retrieve top PubMed studies.",
      class="subtitle")
  )),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      fileInput("gene_file", "Upload Gene List (.txt or .csv)",
                accept=c("text/plain","text/csv",".txt",".csv")),
      helpText("One gene per row; gene names in the first column."),
      
      hr(),
      h5("Search Options"),
      selectInput("sort_by", "Sort results by",
                  choices=c("Relevance"="relevance","Most Recent"="most_recent"),
                  selected="relevance"),
      numericInput("top_n", "Top N results per gene per keyword",
                   value=3, min=1, max=10, step=1),
      checkboxInput("include_total", "Add 'Total Hits' column", value=TRUE),
      
      hr(),
      h5("Date Filter (optional)"),
      fluidRow(
        column(6, numericInput("year_from", "From year", value=NA,
                               min=1900, max=2100, step=1)),
        column(6, numericInput("year_to",   "To year",   value=NA,
                               min=1900, max=2100, step=1))
      ),
      helpText("Leave blank for no date restriction."),
      
      hr(),
      actionButton("add_combo", "+ Add keyword set",
                   class="btn-primary btn-block"),
      helpText("Each box: [Gene] AND [your keywords]."),
      br(),
      uiOutput("keyword_ui"),
      
      # Live debug: shows what keywords are currently active
      div(class="debug-box", textOutput("active_kw_debug")),
      
      hr(),
      numericInput("delay", "Delay between queries (seconds)",
                   value=0.35, min=0.1, max=5, step=0.05),
      helpText("Keep ≥ 0.34 s to respect NCBI rate limits."),
      
      hr(),
      actionButton("run_btn", "Run PubMed Query",
                   icon=icon("search"), class="btn-primary btn-block"),
      br(),
      div(id="status_box", textOutput("status_text")),
      
      hr(),
      downloadButton("download_btn",      "Download CSV",          class="btn-block"),
      br(),
      downloadButton("download_xlsx_btn", "Download Excel (.xlsx)", class="btn-block")
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
  
  # Track keyword slots: a list of named entries so order is explicit
  # Each entry: list(slot_id = "kw_1", label = "Keyword set 1")
  kw_slots     <- reactiveVal(list(
    list(slot_id="kw_1", label="Keyword Set 1")
  ))
  next_slot_id <- reactiveVal(2L)
  
  # Add a new keyword slot
  observeEvent(input$add_combo, {
    sid    <- paste0("kw_", next_slot_id())
    label  <- paste("Keyword Set", next_slot_id())
    slots  <- kw_slots()
    kw_slots(c(slots, list(list(slot_id=sid, label=label))))
    next_slot_id(next_slot_id() + 1L)
  })
  
  # Remove a slot
  observeEvent(input$remove_kw, {
    sid   <- input$remove_kw
    slots <- kw_slots()
    kw_slots(Filter(function(s) s$slot_id != sid, slots))
  })
  
  # Render keyword input boxes
  output$keyword_ui <- renderUI({
    slots <- kw_slots()
    tagList(lapply(slots, function(s) {
      sid <- s$slot_id
      div(class="kw-block",
          tags$button("✕", class="remove-btn",
                      onclick=sprintf(
                        "Shiny.setInputValue('remove_kw','%s',{priority:'event'})", sid
                      )
          ),
          div(class="kw-label", s$label),
          textInput(
            inputId     = paste0("text_", sid),
            label       = NULL,
            placeholder = "e.g. ovarian AND cancer AND therapy"
          )
      )
    }))
  })
  
  # Read active keywords — keyed list: slot_id -> keyword text
  # Isolated so we only capture at run-time, not reactively mid-type
  read_active_keywords <- function() {
    slots <- kw_slots()
    result <- list()
    for (s in slots) {
      val <- trimws(input[[paste0("text_", s$slot_id)]] %||% "")
      if (nchar(val) > 0) {
        result[[s$slot_id]] <- list(text=val, label=s$label)
      }
    }
    result
  }
  
  # Debug display so user can see what's registered
  output$active_kw_debug <- renderText({
    # Reactive dependency on all text_ inputs
    slots <- kw_slots()
    vals  <- vapply(slots, function(s) {
      v <- trimws(input[[paste0("text_", s$slot_id)]] %||% "")
      if (nchar(v) > 0) paste0(s$label, ": ", substr(v,1,30))
      else paste0(s$label, ": (empty)")
    }, character(1))
    paste("Active:", paste(vals, collapse=" | "))
  })
  
  # Gene list
  gene_list <- reactive({
    req(input$gene_file)
    ext <- tolower(tools::file_ext(input$gene_file$name))
    df  <- if (ext == "csv") {
      read.csv(input$gene_file$datapath, header=TRUE, stringsAsFactors=FALSE)
    } else {
      tryCatch(
        read.table(input$gene_file$datapath, header=TRUE,
                   stringsAsFactors=FALSE, sep="\t"),
        error = function(e)
          data.frame(Gene=readLines(input$gene_file$datapath),
                     stringsAsFactors=FALSE)
      )
    }
    genes <- trimws(as.character(df[[1]]))
    genes[nchar(genes) > 0]
  })
  
  output$gene_preview <- renderPrint({
    genes <- gene_list()
    cat(sprintf("Total genes loaded: %d\n\n", length(genes)))
    cat(paste(head(genes, 30), collapse="\n"))
    if (length(genes) > 30)
      cat(sprintf("\n... and %d more", length(genes) - 30))
  })
  
  date_filter <- reactive({
    yf <- input$year_from
    yt <- input$year_to
    if (!is.na(yf) && !is.na(yt)) sprintf(" AND %d:%d[pdat]", yf, yt)
    else if (!is.na(yf))          sprintf(" AND %d:3000[pdat]", yf)
    else if (!is.na(yt))          sprintf(" AND 1900:%d[pdat]", yt)
    else ""
  })
  
  # Snapshot of keywords at run time — stored so summary tab can use it
  run_keywords <- reactiveVal(NULL)
  
  results <- eventReactive(input$run_btn, {
    genes <- gene_list()
    kws   <- read_active_keywords()   # list: slot_id -> list(text, label)
    topn  <- isolate(input$top_n)
    delay <- isolate(input$delay)
    df    <- isolate(date_filter())
    srt   <- isolate(input$sort_by)
    inc_total <- isolate(input$include_total)
    
    validate(
      need(length(genes) > 0, "Please upload a gene list first."),
      need(length(kws)   > 0, "Please enter at least one keyword set.")
    )
    
    run_keywords(kws)   # save for summary tab
    
    out <- data.frame(Gene=genes, stringsAsFactors=FALSE)
    
    kw_ids <- names(kws)
    
    for (j in seq_along(kw_ids)) {
      sid      <- kw_ids[j]
      kw_text  <- kws[[sid]]$text
      kw_label <- kws[[sid]]$label
      # Use the actual keyword text as the column prefix, sanitised for column names
      kw_clean   <- gsub("[^A-Za-z0-9]+", "_", kw_text)
      kw_clean   <- gsub("^_|_$", "", kw_clean)
      kw_clean   <- substr(kw_clean, 1, 40)          # cap length
      col_prefix <- paste0(kw_clean)
      
      top_cols <- paste0(col_prefix, "_Top", seq_len(topn))
      hit_col  <- paste0(col_prefix, "_TotalHits")
      
      out[top_cols] <- NA_character_
      if (inc_total) out[[hit_col]] <- NA_integer_
      
      message(sprintf("[Query] Set %d (%s): '%s'", j, sid, kw_text))
      
      withProgress(
        message = sprintf("Keyword set %d/%d: %s", j, length(kws),
                          substr(kw_text, 1, 30)),
        value = 0, {
          
          for (i in seq_along(genes)) {
            q   <- paste0(genes[i], " AND (", kw_text, ")", df)
            top <- get_top_pubmed_results(q, n=topn, sort_by=srt)
            
            message(sprintf("  Gene=%s  hits=%d  total=%d",
                            genes[i], length(top$ids), top$total_count))
            
            for (k in seq_len(topn)) {
              if (length(top$links) >= k && !is.na(top$links[k]))
                out[i, top_cols[k]] <- top$links[k]
            }
            
            if (inc_total) out[i, hit_col] <- top$total_count
            
            incProgress(1/length(genes),
                        detail=sprintf("%d/%d — %s", i, length(genes), genes[i]))
            Sys.sleep(delay)
          }
        }
      )
    }
    out
  })
  
  output$status_text <- renderText({
    df  <- results()
    kws <- run_keywords()
    n_populated <- sum(!is.na(df[, grep("_Top", names(df)), drop=FALSE]))
    sprintf("Done — %d genes × %d keyword set(s) | %d links found.",
            nrow(df), length(kws), n_populated)
  })
  
  make_prefix <- function(kw_text) {
    p <- gsub("[^A-Za-z0-9]+", "_", kw_text)
    p <- gsub("^_|_$", "", p)
    substr(p, 1, 40)
  }
  
  output$col_visibility_ui <- renderUI({
    kws <- run_keywords()
    if (is.null(kws) || length(kws) <= 1) return(NULL)
    pfxs   <- vapply(kws, function(k) make_prefix(k$text), character(1))
    labels <- paste0("Set ", seq_along(kws), ": ", substr(vapply(kws, `[[`, character(1), "text"), 1, 35))
    checkboxGroupInput("visible_searches", "Show keyword sets:",
                       choices  = setNames(pfxs, labels),
                       selected = pfxs,
                       inline   = TRUE)
  })
  
  filtered_results <- reactive({
    df  <- results()
    vis <- input$visible_searches
    if (is.null(vis)) return(df)
    pattern <- paste(paste0("^", vis, "_"), collapse="|")
    keep    <- c("Gene", grep(pattern, names(df), value=TRUE))
    df[, intersect(keep, names(df)), drop=FALSE]
  })
  
  output$results_table <- renderDT({
    df <- filtered_results()
    top_idx <- grep("_Top", names(df)) - 1L
    datatable(
      df,
      rownames = FALSE,
      escape   = FALSE,
      filter   = "top",
      options  = list(
        pageLength = 25,
        scrollX    = TRUE,
        columnDefs = if (length(top_idx))
          list(list(width="280px", targets=top_idx))
        else list()
      )
    )
  })
  
  # Summary tab
  summary_df <- reactive({
    df  <- results()
    kws <- run_keywords()
    req(!is.null(kws))
    
    # Rebuild the same kw_clean prefix used during query
    make_prefix <- function(kw_text) {
      p <- gsub("[^A-Za-z0-9]+", "_", kw_text)
      p <- gsub("^_|_$", "", p)
      substr(p, 1, 40)
    }
    
    data.frame(
      Keyword_Set = paste0("Set ", seq_along(kws)),
      Keywords    = vapply(kws, `[[`, character(1), "text"),
      Genes_With_Any_Hit = vapply(seq_along(kws), function(j) {
        pfx  <- make_prefix(kws[[j]]$text)
        cols <- grep(paste0("^", pfx, "_Top"), names(df), value=TRUE)
        if (!length(cols)) return(0L)
        sum(rowSums(!is.na(df[cols, drop=FALSE])) > 0)
      }, integer(1)),
      Total_Links = vapply(seq_along(kws), function(j) {
        pfx  <- make_prefix(kws[[j]]$text)
        cols <- grep(paste0("^", pfx, "_Top"), names(df), value=TRUE)
        if (!length(cols)) return(0L)
        sum(!is.na(df[cols, drop=FALSE]))
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
    barplot(
      sm$Genes_With_Any_Hit,
      names.arg = paste0(sm$Keyword_Set, "\n(",
                         substr(sm$Keywords, 1, 30), ")"),
      horiz=TRUE, las=1, col="#3498db", border=NA,
      main="Genes with ≥1 PubMed hit per keyword set",
      xlab="Number of genes"
    )
  })
  
  # ── Download helpers ─────────────────────────────────────────────────────
  # Extract title text and URL from stored HTML anchor tags
  html_title <- function(x)
    gsub("<a href='[^']+' target='_blank'[^>]*>([^<]+)</a>", "\\1", x, perl=TRUE)
  html_url <- function(x)
    gsub("<a href='([^']+)'[^>]*>.*?</a>", "\\1", x, perl=TRUE)
  
  # For CSV: replace each _TopN HTML column in-place with just the title text.
  # CSV can't embed clickable links, so the title is plain text only.
  prepare_csv <- function(df) {
    link_cols <- grep("_Top\\d+$", names(df), value=TRUE)
    for (col in link_cols) {
      titles     <- html_title(df[[col]])
      df[[col]]  <- ifelse(is.na(df[[col]]), NA_character_, titles)
    }
    df
  }
  
  # For Excel: replace each _TopN column in-place with =HYPERLINK(url,"Title")
  # so there is exactly one column per result, showing the title as a blue link.
  prepare_xlsx <- function(df) {
    link_cols <- grep("_Top\\d+$", names(df), value=TRUE)
    # Store titles/urls keyed by col name before we overwrite anything
    meta <- lapply(link_cols, function(col) {
      list(
        titles = html_title(df[[col]]),
        urls   = html_url(df[[col]])
      )
    })
    names(meta) <- link_cols
    list(df=df, meta=meta, link_cols=link_cols)
  }
  
  output$download_btn <- downloadHandler(
    filename = function() paste0("pubmed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file)
      write.csv(prepare_csv(results()), file, row.names=FALSE, na="")
  )
  
  output$download_xlsx_btn <- downloadHandler(
    filename = function() paste0("pubmed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"),
    content  = function(file) {
      if (!requireNamespace("openxlsx", quietly=TRUE)) {
        showNotification("Install openxlsx: install.packages('openxlsx')", type="error")
        return()
      }
      library(openxlsx)
      
      prep <- prepare_xlsx(results())
      df   <- prep$df
      meta <- prep$meta
      
      # Replace HTML in link columns with plain title text for writeDataTable,
      # then overwrite each cell with =HYPERLINK() formula afterwards.
      for (col in prep$link_cols) {
        df[[col]] <- meta[[col]]$titles   # plain text placeholder
      }
      
      wb <- createWorkbook()
      addWorksheet(wb, "PubMed Results")
      writeDataTable(wb, 1, df, tableStyle="TableStyleMedium9")
      setColWidths(wb, 1, cols=seq_len(ncol(df)), widths="auto")
      
      link_style <- createStyle(fontColour="#2980b9", textDecoration="underline")
      
      for (col in prep$link_cols) {
        col_idx <- which(names(df) == col)
        titles  <- meta[[col]]$titles
        urls    <- meta[[col]]$urls
        for (row in seq_len(nrow(df))) {
          url   <- urls[row];   if (is.na(url)   || nchar(url)   == 0) next
          title <- titles[row]; if (is.na(title) || nchar(title) == 0) next
          title_esc <- gsub('"', '""', title)
          writeFormula(wb, 1,
                       x        = sprintf('HYPERLINK("%s","%s")', url, title_esc),
                       startRow = row + 1, startCol = col_idx)
          addStyle(wb, 1, style=link_style, rows=row+1, cols=col_idx)
        }
      }
      saveWorkbook(wb, file, overwrite=TRUE)
    }
  )
}

shinyApp(ui=ui, server=server)