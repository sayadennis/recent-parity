use NM_BI;

DECLARE @cohort_id AS INT = 5096;

IF object_id('tempdb..#COHORTDEMO') IS NOT NULL
	DROP TABLE #COHORTDEMO;

IF object_id('tempdb..#BREASTCANCER') IS NOT NULL
	DROP TABLE #BREASTCANCER;

SELECT DISTINCT p.patient_key, p.ir_id, p.west_mrn AS EPIC_mrn, p.nmh_mrn AS Powerchart_mrn, p.first_name, p.last_name, p.birth_date, YEAR(de.start_datetime) AS diagnosis_year
INTO #BREASTCANCER
FROM dim.vw_patient_current p
	JOIN FSM_Analytics.[fsm_ids_dm].[patient_consent_status] con
		ON p.ir_id = con.ir_id
		AND con.signed_consent = 1
		AND con.revoked_consent = 0
	INNER JOIN fact.vw_diagnosis_event de
		ON p.patient_key = de.patient_key
		AND  p.is_test_patient = 0
		AND p.full_name NOT LIKE '%test%'
		AND p.age >=18
		AND p.gender = 'F' -- Make sure they are women
		AND de.source_system_key IN (1,2,3) -- clarity, NMFF_Clarity, Cerner
		AND YEAR(de.start_datetime) >= 2010
		AND YEAR(de.start_datetime) <= 2020
		AND DATEDIFF(YY, de.start_datetime, p.birth_date) <= 50
	INNER JOIN dim.vw_diagnosis_event_profile dp
		ON de.diagnosis_event_profile_key = dp.diagnosis_event_profile_key
		AND dp.event_type = 'Encounter Diagnosis'	-- For Cancer
	INNER JOIN dim.vw_diagnosis_terminology term
		ON de.diagnosis_key = term.diagnosis_key
		AND term.diagnosis_code_set IN ('ICD-9-CM', 'ICD-10-CM')
		AND (term.diagnosis_code IN  ('174.0', '174.1', '174.2', '174.3', '174.4', '174.5', '174.8', '174.9', '233.0')
			OR term.diagnosis_code_base IN ('C50', 'D05')
			);


SELECT DISTINCT p.ir_id, p.west_mrn AS EPIC_mrn, p.nmh_mrn AS Powerchart_mrn, p.first_name, p.last_name, p.birth_date, bc.diagnosis_year
INTO #COHORTDEMO
FROM dim.vw_patient_current p
	JOIN fsm_analytics.cohort.cohort_patients cp 
		ON cp.source_ir_id = p.ir_id 
            AND cp.is_dltd_flg = 0 
            AND cp.cohort_id = @cohort_id 
	JOIN fsm_analytics.cohort.cohorts c 
        ON c.cohort_id = cp.cohort_id 
            AND c.is_dltd_flg = 0
	JOIN #BREASTCANCER bc
		ON bc.patient_key = p.patient_key
;

SELECT DISTINCT *
FROM #COHORTDEMO;

