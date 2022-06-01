# to run sometimes

# to tidy the description
usethis::use_tidy_description()

# to rebuild the documentation
devtools::document()

# minor performance penalty using :: on the order of 5µs, so it will only matter if you call the function millions of times
