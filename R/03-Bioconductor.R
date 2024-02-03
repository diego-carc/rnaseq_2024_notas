# ### Bioconductor
#
# Bioconductor is an open-source project for the development and sharing of R bioinformatic tools useful in the analysis of massive biological data.
#
# The core team in Bioconductor has the task of developing and maintaining packages and optimizing the use of computational resources.
#
# References.
#
# ### How to find packages
#
# There are several kinds of packages available in Bioconductor:
#
#   - Software: Statistics tools
# - Annotation: Curated information obtained from databases
# - Experiment: Share experimental data
# - Workflows: Documentation packages describing the steps to make a certain type of analysis and involve several other packages.
# - Books: Meta-documentation like describing how to make an analysis using packages available in Bioconductor.
#
# BiocViews is a tree structure organizing the packages available in Bioconductor. We can find packages related to a certain type of study or we can search directly by authors.
#
# ### Bioconductor packages structure
#
# The website for a Bioc package contains a few tags that show us the platforms available to download the package, some statistics about the usage, and the state of maintenance.
#
#The **Installation** section contains a code bloc with the command to download the package

# ```r
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("recount3")
# ```

# **Documentation** contains at least a file (.pdf or .html) describing the proper usage of the package.
#
# Bioconductor owns a git server. Bioc Package Browser allows us to check the package source code.
#
# ### Bioc Branches
#
# **Release:** The stable version of a package
#
# **Devel:** The experimental version of a package
#
# Bioconductor is updated every 6 months.
#
# The community of Bioconductor is currently active in Conferences and resources.

### Exercise
# We discussed the packages "GenomicPlot" and "MethPed".
# We noticed that the newest package is "GenomicPlot" with less than 6 months in
# Bioconductor.This package has verification for all platforms and its build state is
# "ok". Its rank is modest but it is likely due to its short life. Also, the last update
# was during the last week, so the package seems to be currently maintained.
# On the other hand, "MethPad" is a 7.5-year-old package ranked 1797/2266. It is
# verified for all platforms and has a build state  "ok". It has no opened questions in
# the support section, so we argued whether this means the package has no troubles or the
# package is not being used enough. We concluded that the package seems to be robust
# enough for proper maintenance.
#
