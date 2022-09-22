# Manually adding patients to our cohort via extensive chart-review

## Background 

Germline genetic testing results are first reported on a PDF file sent from the vendor. Typically, the genetic counselor then enters this information onto the EHR systm when it is received. However, it appeared that this was not always the case - i.e. some patients appeared highly likely to have genetic testing results available (their notes would end with "genetic testing results pending") although this information was not available in the database. 

## Getting around this issue 

With the extensive help of our physician collaborators, genetic testing and other data of these "patients missing genetic testing data" was obtained through manual chart-review and data entry. 

## Scripts descriptions 

Scripts 01-04 describes the process prior to chart review, obtaining as much gynecological and pathological information from the EHR for this patient set as possible prior to chart review. Script 05 joins the reviewed data with the previously obtained dataset and the following scripts perform some preliminary analysis with the combined data. 
