# Resubmission

This is a resubmission. In this version I have:

-   Update Description field in DESCRIPTION file to meet the CRAN requirements.

## R CMD check results

0 errors \| 0 warnings \| 1 note

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

-   We saw 0 new problems
-   We failed to check 0 packages

# PREVIOUS CRAN submissions

## Resubmission for MaxWiK_1.0.4.tar.gz

##### Requests from CRAN team

1.  Please do not start the description with "This package", package name, title or similar.

2.  \\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \\dontrun{} adds the comment ("\# Not run:") as a warning for the user. Does not seem necessary. Please replace \\dontrun with \\donttest. Please unwrap the examples if they are executable in \< 5 sec, or replace \\dontrun{} with \\donttest{}.

3.  Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies. Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir(). -\> R/utils.R

4.  Please do not modify the global environment (e.g. by using \<\<-) in your functions. This is not allowed by the CRAN policies. -\> R/Sampler_iKernel.R

##### Answers on requests

1.  Description is changed with recommended format like 'MaxWiK' is a machine learning method beased on maxima weighted Isolation Kernel approach, etc.

2.  \\dontrun{} is deleted as not needed in all the examples. Instead tempdir() is used as input folder.

3.  The default path is changed to NULL in the function MaxWiK_templates(). The examples in the vignettes are modified to use tempdir() function.

4.  The code of the function Sampler_iKernel() is modified to avoid using the global variables. Previously it was a counter of the number of iterations warapped in the function inside the other function. Now it is a parameter of the function.

## Resubmission for MaxWiK_1.0.3.tar.gz

This is a resubmission. The package MaxWiK_1.0.2.tar.gz does not pass the incoming checks automatically because of two notes.

This version fixed the notes:

##### Automatical check results

2 NOTES:

1.  

-   checking CRAN incoming feasibility ... [4s/6s] NOTE Maintainer: ‘Yuri Nagornov [nagornov.yuri\@gmail.com](mailto:nagornov.yuri@gmail.com){.email}’

New submission

Possibly misspelled words in DESCRIPTION: Iurii (6:667, 6:728) Nagornov (6:673, 6:734)

2.  

-   checking for non-standard things in the check directory ... NOTE Found the following files/directories: ‘MaxWiK.ABC.R’ ‘MaxWiK.Predictor.R’ ‘MaxWiK.Sampling.R’

##### Answers on NOTES:

1.  

Iurii and Nagornov are the correct spelling of the author's name. Actually I use Yuri and Iurii spelling in different cases. To fix it I have updated inst/WORDLIST file.

And this is the submission of a new package.

2.  

Files ‘MaxWiK.ABC.R’ ‘MaxWiK.Predictor.R’ ‘MaxWiK.Sampling.R’ are used as templates for demonstration of the machine learning method for three different problems. They are described in the vignette. To avoid copying of these files to the user's working directory, I have added "\\dontrun{}" to the manual, and also modified 'MaxWiK_templates( dir = tempdir() )' examples in the vignettes using tempdir() function.
