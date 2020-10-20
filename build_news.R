build_news <- function(outfile) {

  pavo_news <- readLines("https://raw.githubusercontent.com/rmaia/pavo/master/NEWS.md")
  lightr_news <- readLines("https://raw.githubusercontent.com/ropensci/lightr/main/NEWS.md")

  # We only display the changes for the two latest versions
  pavo_news <- pavo_news[seq_len(grep("^# ", pavo_news)[3]-1)]
  lightr_news <- lightr_news[seq_len(grep("^# ", lightr_news)[3]-1)]

  # Decrease headings level by one
  pavo_news <- gsub("^(#+)", "\\1#", pavo_news)
  lightr_news <- gsub("^(#+)", "\\1#", lightr_news)

  content <- paste(
    "---",
    "layout: page",
    "title: Updates",
    "---",
    "",
    "The most recent package updates:",
    "",
    paste(pavo_news, collapse = "\n"),
    "---",
    paste(lightr_news, collapse = "\n"),
    sep = "\n"
  )

  write(content, outfile)

}
