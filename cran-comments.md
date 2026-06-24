## CRAN comments for prova 0.7.0

First submission.

### R CMD check results

Linux (R-devel), Windows (win-builder, R-devel):

0 errors | 1 warnings | 3 notes

Please see explanations below about the Warning and Notes.


### Reverse dependency check:

No reverse dependencies.


### Explanation about the Warning

The Warning appears because of the way the 'Nimble' package is handled, and is motivated as follows. I'd be grateful if the CRAN maintainers could provide an alternative solution, or grant the Warning to be acceptable.

- The functions of this package can become very memory expensive, so the code is written as to minimize memory use.

- The function 'learn()' is the only one that uses the 'Nimble' package. But it only loads that package in one or more separate, non-interactive R sessions, by means of the 'parallel' package. Those sessions are closed before 'learn()' returns its output. This function does not need 'Nimble' at all outside those R sessions.

- The 'Nimble' package is not used or needed by any other functions of the present package.

It is for the reasons above that I would prefer if the present package did not load the package 'Nimble' (with its numerous dependencies) by default. This would use more than 100 MB of memory for no reason. Of course users can still load the 'Nimble' package by themselves if they so wish.

To achieve this situation, 'Nimble' is listed in the "Imports:" field of the DESCRIPTION file of the present package, but it is not called in any 'import()' or 'importFrom()' in the NAMESPACE.

As an alternative solution I considered reporting 'Nimble' only in the 'Suggests:' field. However, that would somewhat be a misrepresentation of the dependency.

Again, I'd be grateful if the CRAN maintainers could provide an alternative solution, or grant the Warning and the handling of the 'Nimble' package to be acceptable.

### Explanation about the Notes

- Note 1: Maintainer has another package on CRAN ('Pinference').

- Note 2: No visible bindings for several global variables (in functions 'learn()' and 'workerfun()' which handles package 'Nimble' in a separate R session): this is because of how some variables are internally handled by the 'Nimble' package.

- Note 3: (a) Example for function 'learn()' takes around 30 s to 60 s. Unfortunately it is not possible to further reduce the computation time of this function. If necessary I can put the example in a '\donttest{}' context. (b) Examples for function 'qPr()' take around 10 s to 30 s. I can reduce the number of examples if this is considered too much.

