use NM_BI;

/***
The purpose of this code is to confirm the breast cancer diagnoses that we want to include in our cohort.
***/

SELECT DISTINCT term.diagnosis_code, term.diagnosis_name
FROM dim.vw_diagnosis_terminology term
WHERE term.diagnosis_code_set IN ('ICD-9-CM', 'ICD-10-CM')
		AND (term.diagnosis_code IN  ('174.0', '174.1', '174.2', '174.3', '174.4', '174.5', '174.8', '174.9', '233.0')
			OR term.diagnosis_code IN ('C50', 'D05')
			OR term.diagnosis_code LIKE 'C50.[0-9]'
			OR term.diagnosis_code LIKE 'C50.[0-9]1%'
			OR term.diagnosis_code LIKE 'D05.[1-9]%'
		);
