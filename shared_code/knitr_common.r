knitr.project_name <- "Early embryo Pol II"

figure_path <- function(filename="") {
  paste(knitr.figure_dir, "/", filename, sep="")
}

# Output HTML table
html_table <- function(df, row.names=FALSE, col.names=TRUE) {
  print(xtable(df), type="html", include.rownames=row.names, include.colnames=col.names)
}

# Format number with commas
pn <- function(i) {
  prettyNum(i, big.mark=",")
}

# Force knitr to stop evaluation when an error is encountered
opts_knit$set(stop_on_error=2L)

# Wrap output and code
options(width=80)

# Don't echo code
opts_chunk$set(echo=FALSE)

# Don't reformat code
opts_chunk$set(tidy=FALSE)

# Set up figure defaults 
opts_chunk$set(fig.width=7, fig.height=5, fig.path=figure_path())

if(!file.exists(knitr.figure_dir)) dir.create(knitr.figure_dir)
