## CRAN comments for prova 0.7.0

Second submission.

- Changed dependence type from 'Nimble' package to avoid Warning.
- Solved "no visible bindings for several global variables" Note.
- Put examples that take too long computation within '\donttest{}'.


### R CMD check results

Linux (R-devel), Windows (win-builder, R-devel):

0 errors | 0 warnings | 2 notes

Please see explanations below about Notes.


### Reverse dependency check:

No reverse dependencies.


### Explanation about Notes and other items

- Words spellings have been double-checked and are correct.

- URLs to 'dictionary.cambridge.org' are valid; probably they are blocked by Cloudflare when checked by a bot.


- Note 1: Maintainer has another package on CRAN ('Pinference').


- Note 2: The package 'nimble', used by the 'learn()' function of the present package, cannot be only loaded with 'requireNamespace()' and called with '::', as opposed to attached with 'require()' or 'library()'. This is because of the way it works; see discussion at <https://groups.google.com/g/nimble-users/c/dDNE3L_sPxI/m/IyRgmXOaBwAJ>.

However, the function 'learn()' is the only one that uses the 'nimble' package, and *it only uses that package in one or more temporary, separate, non-interactive R sessions*, spawned by the 'parallel' package. Those sessions are closed before 'learn()' returns its output. This function does not need 'nimble' at all outside those temporary R sessions. And the 'nimble' package is *not* used or needed by any other functions of the present package.

The functions of this package can become very memory expensive, so the code is written as to minimize memory use. If 'nimble' were added among the "Imports", it would always be attached upon attaching the present package, unnecessarily occupying around 100 MB of memory.

For this reason I would prefer to leave 'nimble' among the "Suggested" packages, rather than among the "Imports" ones, and attach it with 'require()' the temporary, non-interactive R session where it is used.

I am happy to hear about better solutions from the CRAN reviewers.
