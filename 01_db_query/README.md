# Step 1: Database Querying 

This folder stores the code for querying the necessary data using SSMS. 

0. Determine the method of cohort identification by testing various approaches (`explore_cohort`) and assessing the sensitivity/specificity in comparison to a maually-curated patient list. 
1. Export the internal IDs of the selected cohort. (`cohort_selection.sql`)
2. Query the data necessary to obtain the demographic (`demographic.sql`), gynecological (`notes_initial_counsel.sql`), germline genetic-testing (`notes_genetic_counsel.sql`), and pathological (`notes_pathology.sql`) data of this cohort. 

Note: `diagnosis_codes_names.sql` is used for a sanity-check to confirm the diagnoses we are including. 