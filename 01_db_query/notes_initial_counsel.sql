use NM_BI;

IF object_id('tempdb..#COUNSELNOTES') IS NOT NULL
	DROP TABLE #COUNSELNOTES;

SELECT n.* --n.note_text  
INTO #COUNSELNOTES
FROM fact.vw_note n
	INNER JOIN dim.vw_department dep
		ON n.department_key = dep.department_key
		AND dep.department_specialty = 'Breast Surgery'
		AND n.source_system_table = 'hno_note_text'
WHERE (n.note_text LIKE '%initial consult%'
	OR n.note_text LIKE '%consultation%'
	OR n.note_text LIKE '%new patient history and physical%');

--IF object_id('tempdb..#GENETICSNOTES') IS NOT NULL
--	DROP TABLE #GENETICSNOTES;

--SELECT n.* --n.note_text  
--INTO #GENETICSNOTES
--FROM fact.vw_note n
--	INNER JOIN dim.vw_department dep
--		ON n.department_key = dep.department_key
--		AND (
--			dep.department_specialty = 'Genetics' 
--			OR 
--			dep.department_specialty = 'Genetic Counseling'
--		)
--WHERE n.note_text LIKE '%genetic result addendum%';

--IF object_id('tempdb..#BREASTCANCER') IS NOT NULL
--	DROP TABLE #BREASTCANCER;

--SELECT DISTINCT p.patient_key, p.ir_id, p.west_mrn AS EPIC_mrn, p.nmh_mrn AS Powerchart_mrn, p.first_name, p.last_name, p.birth_date
--INTO #BREASTCANCER
--FROM dim.vw_patient_current p
--	JOIN FSM_Analytics.[fsm_ids_dm].[patient_consent_status] con
--		ON p.ir_id = con.ir_id
--		AND con.signed_consent = 1
--		AND con.revoked_consent = 0
--	INNER JOIN fact.vw_diagnosis_event de
--		ON p.patient_key = de.patient_key
--		AND  p.is_test_patient = 0
--		AND p.full_name NOT LIKE '%test%'
--		AND p.age >=18
--		AND p.gender = 'F' -- Make sure they are women
--		AND de.source_system_key IN (1,2,3) -- clarity, NMFF_Clarity, Cerner
--		AND YEAR(de.start_datetime) >= 2010
--		AND YEAR(de.start_datetime) <= 2020
--		AND DATEDIFF(YY, de.start_datetime, p.birth_date) <= 50
--	INNER JOIN dim.vw_diagnosis_event_profile dp
--		ON de.diagnosis_event_profile_key = dp.diagnosis_event_profile_key
--		AND dp.event_type = 'Encounter Diagnosis'	-- For Cancer
--	INNER JOIN dim.vw_diagnosis_terminology term
--		ON de.diagnosis_key = term.diagnosis_key
--		AND term.diagnosis_code_set IN ('ICD-9-CM', 'ICD-10-CM')
--		AND (term.diagnosis_code IN  ('174.0', '174.1', '174.2', '174.3', '174.4', '174.5', '174.8', '174.9', '233.0')
--			OR term.diagnosis_code_base IN ('C50', 'D05')
--			);

--SELECT DISTINCT bc.ir_id, bc.EPIC_mrn, bc.Powerchart_mrn, bc.first_name, bc.last_name, bc.birth_date, 
--	cn.note_key, cn.source_system_id, cn.source_system_key, cn.source_system_table, cn.created_datetime, cn.updated_datetime --, cn.note_text
--FROM #BREASTCANCER bc
--	JOIN #COUNSELNOTES cn
--		ON cn.patient_key = bc.patient_key
--	JOIN #GENETICSNOTES gn
--		ON gn.patient_key = bc.patient_key;


DECLARE @cohort_id AS INT = 5096;

IF object_id('tempdb..#COHORTCOUNSELNOTES') IS NOT NULL
	DROP TABLE #COHORTCOUNSELNOTES;

SELECT DISTINCT p.ir_id, p.west_mrn AS EPIC_mrn, p.nmh_mrn AS Powerchart_mrn, p.first_name, p.last_name, p.birth_date,
		cn.note_key, cn.source_system_id, cn.source_system_key, cn.source_system_table, YEAR(cn.created_datetime) AS created_year, YEAR(cn.updated_datetime) as updated_year, cn.note_text
INTO #COHORTCOUNSELNOTES
FROM dim.vw_patient_current p
	JOIN fsm_analytics.cohort.cohort_patients cp 
		ON cp.source_ir_id = p.ir_id 
            AND cp.is_dltd_flg = 0 
            AND cp.cohort_id = @cohort_id 
	JOIN fsm_analytics.cohort.cohorts c 
        ON c.cohort_id = cp.cohort_id 
            AND c.is_dltd_flg = 0
	JOIN #COUNSELNOTES cn
		ON p.patient_key = cn.patient_key
;

SELECT *
FROM #COHORTCOUNSELNOTES;
