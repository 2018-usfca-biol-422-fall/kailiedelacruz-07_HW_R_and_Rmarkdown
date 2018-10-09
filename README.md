# R and Rmarkdown homework assignment for BIOL422/423
## Due as a pull request on Friday, October 16, 2018 before 11:59 pm

The **goal of this assignment** is to introduce you to the use of R, RStudio, and Rmarkdown to analyze data. The input data for this week's analyses will be the BLAST output files from the BLAST scripting homework. 

The data we will be using is from the NCBI Sequence Read Archive study number ERP022657. A summary of the information is available [here](https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_128047291_130.14.22.33_5555_1505945515_1626731749_0MetA0_S_HStore&query_key=5). The metadata from this study is included in the git repository for this assignment in a `data/metadata` directory.

Here's the abstract from the [original study](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=ERP022657) by Fierer et al.:

> Recent work has demonstrated that the diversity of skin-associated bacterial communities is far higher than previously recognized, with a high degree of interindividual variability in the composition of bacterial communities. Given that skin bacterial communities are personalized, we hypothesized that we could use the residual skin bacteria left on objects for forensic identification, matching the bacteria on the object to the skin-associated bacteria of the individual who touched the object. Here we describe a series of studies demonstrating the validity of this approach. We show that skin-associated bacteria can be readily recovered from surfaces (including single computer keys and computer mice) and that the structure of these communities can be used to differentiate objects handled by different individuals, even if those objects have been left untouched for up to 2 weeks at room temperature. Furthermore, we demonstrate that we can use a high-throughput pyrosequencing-based approach to quantitatively compare the bacterial communities on objects and skin to match the object to the individual with a high degree of certainty. Although additional work is needed to further establish the utility of this approach, this series of studies introduces a forensics approach that could eventually be used to independently evaluate results obtained using more traditional forensic practices.

For this assignment, you will only need to work in one file, the `Rmd` file entitled `Analysis_of_BLAST_Results.Rmd`. There is some starter code in that file, with comments on how the various functions work. You should modify this document as needed, and add code chunks and markdown text to complete the following set of tasks.

This assignment very much builds on the work you did last week with the BLAST script and your analysis of the results. In this assignment we are digging deeper into the results using R. Just like last time, you should use your results as a launching off point to do some additional research on the taxa you find. You can use Google Scholar or some other tool to search the peer-reviewed literature. 

This report could focus in on, for example: different bacterial taxa on male vs female hands, or the similarity of mouse surfaces to student hands, or just what is found on computer mice, etc. It matters less to me what specific aspect of the data you focus on, and more important that you spend some time researching the biological implications of the taxa you identified using the BLAST search and how they are different between different sample types or individuals.

Please follow the instructions carefully and read them all before getting started.

This second assignment will be worth 20 points. The grading breakdown will be as follows (notice it is slightly different from past assignments):

* 10 points - Completes all required steps, including the analysis report and associated literature search (as outlined below)
* 3 points - Scripts are appropriately commented and well organized
* 3 points - Appropriate use of git to version control the steps, including adding and committing the appropriate files at the specific steps below, and writing informative and appropriately formatted commit messages
* 4 points - Pull Request passes automated checks

You must submit your work as a Pull Request to the class organization ('2018-usfca-biol-422-fall') on GitHub by 11:59 pm on Tuesday, October 16 for full credit. Your PR should pass a set of automated checks (more detail in the last step, below). We will also be peer reviewing the code after it is submitted, as usual.

Steps:

1. Fork this repository to your own GitHub account.
1. In a new RStudio session on your laptop (no open files, projects, etc), select "File -> New Project..." from the menu. Then click on "Version Control", then "Git", and then paste the URL from your **forked** repository (the one under your GitHub account, that should include something like YourName/07_HW_R_and_Rmarkdown) into the "Repository URL" field. Then choose where you would like the folder to be on your computer and click "Create Project". This will use `git clone` to download the repository from GitHub and get RStudio set up for you to work on it. This should also properly set up the git remote so that you can push and pull from within RStudio. 
1. Run through the existing code in the document manually to make sure each step works and to build the combined dataset that you will use for analysis. What I mean by this is to execute the code chunks by hand instead of knitting the document. Then you can click on the dataset called `joined_BLAST_data_metadata` in the Environment tab in the top right pane of RStudio to see what the whole dataset looks like.
1. You should include in this analysis report all of the parts of a standard scientific paper, including at a minimum: Introduction, Methods, Results, Discussion. I have added top-level headers by using a single `#`; if you want to make additional subsection headers you can use `##` or `###`. Each of these sections should have a least 1-2 paragraphs of markdown text, briefly outlining some background on the data we are using (you can draw on the article we read by Fierer *et al.* as well as the abstract above), the Methods that were used to generate and analyze the data, and then your Results, and the interpretation of those results in a Discussion section. I think it makes the most sense to have all of the code chunks fall under the Results section. 
1. In the results section, using dplyr's `filter()` function, select a subset of the complete data and visualize it using a ggplot histogram. This histogram could be on any of the numeric columns returned from the BLAST run, but should be something that you think is useful to look at. I highly recommend you download and take a look at the [Data Transformation with dplyr](https://raw.githubusercontent.com/rstudio/cheatsheets/master/data-transformation.pdf) and [R markdown](https://raw.githubusercontent.com/rstudio/cheatsheets/master/rmarkdown-2.0.pdf) cheat sheets while you work. I have an example of this done for you.
1. Repeat this process for at least 3 other different subsets of the data (creating at least 4 figures total). You can select whatever subset and whatever column you like for these, but the different things you choose to examine should be coherent. You can make other types of figures if you wish.
1. After the code that generates each of them, write a short bit of markdown text, right in the results section, interpreting what you see. What is expected? What is surprising? What might be worth following up on next?
1. You can also include and discuss the output of the summary table with counts of sequences for each grouping of subject and sample type (i.e. the `kable(table())` function near the end of the starter Rmd), but it isn't necessary. If you want to discuss it, please include some markdown text doing so, if not, remove the table altogether.
1. Commit the script as you work on it, whenever you make a good chunk of progress. Make sure you write an [appropriate commit message](https://chris.beams.io/posts/git-commit/).
1. After you have finished the analysis, your `Rmd` knits successfully, and you have fixed all errors or warnings (except the 'there is no variable with this name in scope' warnings you see when using dplyr which is still a bug that hasn't been fixed by RStudio), be sure to add a commit marking this milestone (and push back up to GitHub just to be safe!).
1. Once that's all done, add, commit, and push the `Rmd` file itself, **as well as the `md` file and the folder of images that are generated** when you knit back to your fork of the original repository on GitHub using the git interface in RStudio. Remember that you can only push what you have committed, so be sure all of your work is committed. Be sure to save your files often, and check the git tab in RStudio as you work.
2. Submit a Pull Request back to the organization repository to submit your assignment. Make sure the Pull Request (PR) has a useful description of the changes you made. When you submit the pull request, a set of automated tests will run on it, to check that the `Rmd` file knits properly without error and that you have fixed all of the potential code style and syntax errors in your document. If you want to rigorously check for any syntax errors while you are working on the code, you can install the `lintr` package by typing the following at the console:

```r
install.packages("devtools")
devtools::install_github("jimhester/lintr")
```

and then whenever you want to check your code you can run at the R console in RStudio:

```r
lintr::lint(filename = "Analysis_of_BLAST_Results.Rmd")
```

**Pro Tip** Save often, commit often, push often! If you have any questions or need clarification of what it is I'd like you to do, please ask me sooner rather than later so you stay on the right track for completing this assignment on time.

**NOTE:** There are some files in this repository you may not have seen before. They are used for automating testing. 

##### Infrastructure for Automated Software Testing

- `.travis.yml`: a configuration file for automatically running [continuous integration](https://travis-ci.com) checks when you submit your pull request, to verify reproducibility of all `.Rmd` notebooks in the repo.  If all `.Rmd` notebooks can render successfully and pass linting (or code style and syntax checks), then the "Build Status" badge above will be green (`build success`), otherwise it will be red (`build failure`).  
- `.Rbuildignore`: a configuration file telling the build script which files to ignore.
- `DESCRIPTION`: a metadata file for the repository, based on the R package standard. Its main purpose here is as a place to list any additional R packages/libraries needed for any of the `.Rmd` files to run.
- `tests/render_rmds.R`: an R script that is run to execute the above described tests, rendering and linting all `.Rmd` notebooks. 

