# lintr over the workflow scripts. Style only; correctness is the other tiers.
lints <- lintr::lint_dir(
  "workflow/scripts",
  linters = lintr::linters_with_defaults(
    line_length_linter = lintr::line_length_linter(100),
    object_name_linter = NULL,      # cd_ prefixes and .CD are deliberate
    commented_code_linter = NULL,   # code in comments is often the point here
    cyclocomp_linter = NULL
  )
)
if (length(lints)) { print(lints); quit(status = 1) }
cat("lintr: clean\n")
