## Resubmission

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

## R CMD check results

0 errors \| 0 warnings \| 1 note

-   This is a new release.

## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

-   We saw 0 new problems
-   We failed to check 0 packages
