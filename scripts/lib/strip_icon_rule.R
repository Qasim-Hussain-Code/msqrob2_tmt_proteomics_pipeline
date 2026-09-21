## Quarto embeds the complete Bootstrap Icons stylesheet in every HTML
## page. One of its icon classes is named after an AI vendor, and this
## repository is to carry no such name anywhere, so the rule for the
## icon at codepoint f914 is removed from the rendered report. The page
## does not use that icon.
##
## Usage: Rscript scripts/lib/strip_icon_rule.R results/report.html
f <- commandArgs(trailingOnly = TRUE)[1]
x <- readLines(f, warn = FALSE)
pattern <- "\\.bi-[a-z0-9-]+::before \\{ content: \"\\\\f914\"; \\}"
n <- sum(grepl(pattern, x))
x <- gsub(pattern, "", x)
writeLines(x, f)
cat(sprintf("strip_icon_rule: removed %d rule(s) from %s\n", n, basename(f)))
