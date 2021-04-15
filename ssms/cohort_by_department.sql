use NM_BI;

SELECT DISTINCT p.ir_id, p.west_mrn AS EPIC_mrn, p.nmh_mrn AS Powerchart_mrn, p.first_name, p.last_name, p.birth_date
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
			)
;

SELECT DISTINCT bc.*
FROM dim.vw_patient_current p
	INNER JOIN fact.vw_encounter_outpatient eo
		ON eo.patient_key = p.patient_key
	INNER JOIN dim.vw_provider prov
		ON eo.visit_provider_key = prov.provider_key
	INNER JOIN #breastcancer bc
		ON bc.ir_id = p.ir_id
	JOIN dim.vw_department dep
		ON dep.department_key = eo.department_key
WHERE dep.department_name LIKE '%genetics%';
