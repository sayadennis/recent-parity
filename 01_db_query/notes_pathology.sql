use NM_BI;

IF object_id('tempdb..#PATHOLNOTES') IS NOT NULL
	DROP TABLE #PATHOLNOTES;

SELECT n.* --n.note_text
INTO #PATHOLNOTES
FROM fact.vw_note n
	INNER JOIN dim.vw_department dep
		ON n.department_key = dep.department_key
		AND dep.department_specialty = 'Breast Surgery'
WHERE (
	n.note_text LIKE '%surgical pathology final report%'
	OR n.note_text LIKE '%surg path final report%'
	OR n.note_text LIKE '%surgical pathology addendum report%'
	OR n.note_text LIKE '%surg path addendum report%'
	OR n.note_text LIKE '%tumor staging report%'
	OR n.note_text LIKE '%breast cancer tumor markers%'
	OR n.note_text LIKE '%breast cancer staging summary%'
);

IF object_id('tempdb..#GENETICSNOTES') IS NOT NULL
	DROP TABLE #GENETICSNOTES;

SELECT n.* --n.note_text  
INTO #GENETICSNOTES
FROM fact.vw_note n
	INNER JOIN dim.vw_department dep
		ON n.department_key = dep.department_key
		AND (
			dep.department_specialty = 'Genetics' 
			OR 
			dep.department_specialty = 'Genetic Counseling'
		)
WHERE n.note_text LIKE '%genetic result addendum%';

SELECT DISTINCT p.patient_key, p.ir_id, p.west_mrn AS EPIC_mrn, p.nmh_mrn AS Powerchart_mrn, p.first_name, p.last_name, p.birth_date
INTO #BREASTCANCER
FROM dim.vw_patient_current p
	--JOIN FSM_Analytics.[fsm_ids_dm].[patient_consent_status] con
	--	ON p.ir_id = con.ir_id
	--	AND con.signed_consent = 1
	--	AND con.revoked_consent = 0
	INNER JOIN fact.vw_diagnosis_event de
		ON p.patient_key = de.patient_key
		AND  p.is_test_patient = 0
		AND p.full_name NOT LIKE '%test%'
		AND p.age >=18
		AND p.gender = 'F' -- Make sure they are women
		AND de.source_system_key IN (1,2,3) -- clarity, NMFF_Clarity, Cerner
		AND YEAR(de.start_datetime) >= 2010
		AND YEAR(de.start_datetime) <= 2020
		AND (YEAR(de.start_datetime)-YEAR(p.birth_date)) <= 50
	INNER JOIN dim.vw_diagnosis_event_profile dp
		ON de.diagnosis_event_profile_key = dp.diagnosis_event_profile_key
		AND dp.event_type = 'Encounter Diagnosis'	-- For Cancer
	INNER JOIN dim.vw_diagnosis_terminology term
		ON de.diagnosis_key = term.diagnosis_key
		AND term.diagnosis_code_set IN ('ICD-9-CM', 'ICD-10-CM')
		AND (term.diagnosis_code IN  ('174.0', '174.1', '174.2', '174.3', '174.4', '174.5', '174.8', '174.9', '233.0')
			OR term.diagnosis_code IN ('C50', 'D05')
			OR term.diagnosis_code LIKE 'C50.[0-9]'
			OR term.diagnosis_code LIKE 'C50.[0-9]1%'
			OR term.diagnosis_code LIKE 'D05.[1-9]%'
			);

SELECT DISTINCT bc.ir_id, bc.EPIC_mrn, bc.Powerchart_mrn, bc.first_name, bc.last_name, bc.birth_date
FROM #BREASTCANCER bc
	JOIN #PATHOLNOTES pn
		ON pn.patient_key = bc.patient_key
	JOIN #GENETICSNOTES gn
		ON gn.patient_key = bc.patient_key;
