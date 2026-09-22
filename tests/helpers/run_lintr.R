# lintr over the workflow scripts, configured by the repository's .lintr file.
#
# Tuned for correctness, not cosmetics. The linters kept in .lintr catch things
# that change behaviour: an undefined symbol, `&` where `&&` was meant,
# `1:length(x)` on an empty vector, `== NA`. The layout linters are dropped
# because lintr's defaults disagree with this project's style (it wants
# 29-space hanging indents in the standalone blocks, which would make them
# unreadable).
#
# object_usage_linter is off. chromadyn is a script layout, not a package, so
# lib/*.R files legitimately use symbols defined in a sibling file that the
# calling script sources first (plots.R uses theme.R's palettes, for example).
# lintr cannot follow a source() whose path is computed, so it reports every
# one of those as an undefined global. The dependency is enforced at runtime
# instead: lib/plots.R refuses to load if theme.R has not been sourced.
#
# .lintr is DCF, so it cannot carry this explanation itself: a leading comment
# block makes read.dcf() reject the file.
lints <- lintr::lint_dir("workflow/scripts")
if (length(lints)) {
  print(lints)
  cat(sprintf("\n%d lint(s). Fix them, or adjust .lintr if the rule is wrong.\n", length(lints)))
  quit(status = 1)
}
cat("lintr: clean\n")
