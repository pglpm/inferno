# Changelog

## Prova v1.0.0

- Improved documentation, also of internal functions.
- Introduced argument `verbose =` in several function, so as to show
  information to console only if the user so desires. This argument is
  `FALSE` for all functions except
  [`learn()`](https://pglpm.github.io/prova/reference/learn.md).
- Added example of imputation in vignette.
- Rewritten README with description of the package’s main features.
- Cleaned up code.

## Prova v0.9.0

- Improved documentation, also of internal functions.
- Removed unused functions.
- Cleaned up code.

## Prova v0.7.0

- Initial CRAN submission.
- Introduced the [`print()`](https://rdrr.io/r/base/print.html) method
  for probability objects.
- Slightly modified output of
  [`Pr()`](https://pglpm.github.io/prova/reference/Pr.md) and related
  functions.
- Added references.
- Modified vignettes.

## Prova v0.6.5

- Added examples to all functions visible to the user.
- Fixed some bugs (in particular a bug that would allow the input of a
  non-existent nominal value giving erroneous probabilities).
- Corrected typos.

## Prova v0.6.0

Renamed the package from **inferno** to **Prova**. It turns out,
unsurprisingly, that the name “inferno” was already in use for several
inference-related software packages of various kinds.
