use NM_BI;

DECLARE @cohort_id AS INT = 5096;

IF object_id('tempdb..#PATHOLNOTES') IS NOT NULL
	DROP TABLE #PATHOLNOTES;

SELECT DISTINCT p.ir_id, p.patient_key, n.source_system_key, n.source_system_table, 
	n.department_key, dep.department_name, dep.department_external_name, dep.department_specialty, dep.location_name,
	n.created_datetime, n.updated_datetime, n.signed_datetime, n.note_type, n.record_status, n.note_text  
INTO #PATHOLNOTES
FROM dim.vw_patient_current p
	INNER JOIN fsm_analytics.cohort.cohort_patients cp 
		ON cp.source_ir_id = p.ir_id
            AND cp.is_dltd_flg = 0 
            AND cp.cohort_id = @cohort_id 
	INNER JOIN fact.vw_note n
		ON n.patient_key = p.patient_key
		AND YEAR(n.created_datetime)>=2010
		AND YEAR(n.created_datetime)<=2020
	INNER JOIN dim.vw_department dep
		ON n.department_key = dep.department_key
		AND dep.department_specialty IN (
				'Genetics', 'Genetic Counseling', 
				'Obstetrics and Gynecology', 'Obstetrics', 'Gynecology', 'Ob/Gyn',
				'Gynecologic Oncology', 'Gynecology/Oncology', 'Oncology',
				'Hematology/Oncology', 'Hematology and Oncology', 
				'Pathology'
		)
		AND n.source_system_table = 'hno_note_text'
WHERE (
	n.note_text LIKE '%surgical pathology final report%'
	OR n.note_text LIKE '%surg path final report%'
	OR n.note_text LIKE '%surgical pathology addendum report%'
	OR n.note_text LIKE '%surg path addendum report%'
	OR n.note_text LIKE '%tumor staging report%'
	OR n.note_text LIKE '%breast cancer tumor markers%'
	OR n.note_text LIKE '%breast cancer staging summary%'
	OR n.note_text LIKE '%surgpath%'
	OR n.note_text LIKE '%surgical pathology report%'
	OR n.note_text LIKE '%synoptic report%'
	OR n.note_text LIKE '%in-situ breast carcinoma checklist%'
	OR n.note_text LIKE '%surgical pathology report%'
);


select top 200 *
from #PATHOLNOTES;

select department_specialty, note_type, count(*)
from #PATHOLNOTES
group by department_specialty, note_type;

select *
from dim.vw_department dep
where department_name like '%pathol%';


use NM_BI;

DECLARE @cohort_id AS INT = 5096;

IF object_id('tempdb..#TESTPATHOLNOTES') IS NOT NULL
	DROP TABLE #TESTPATHOLNOTES;

SELECT DISTINCT p.ir_id, p.patient_key, p.first_name, p.last_name, p.birth_date, 
	dep.department_name, dep.department_specialty, 
	n.created_datetime, n.updated_datetime, n.signed_datetime, n.note_type, n.record_status, n.note_text  
INTO #TESTPATHOLNOTES
FROM dim.vw_patient_current p
	INNER JOIN fsm_analytics.cohort.cohort_patients cp 
		ON cp.source_ir_id = p.ir_id
            AND cp.is_dltd_flg = 0 
            AND cp.cohort_id = @cohort_id 
	INNER JOIN fact.vw_note n
		ON n.patient_key = p.patient_key
		AND n.source_system_table = 'hno_note_text'
		AND note_type IN ('Progress Notes', 'Letter')
		and n.note_text LIKE '%breast%'
	INNER JOIN dim.vw_department dep
		ON n.department_key = dep.department_key
		AND dep.department_specialty IN ('Genetics', 'Breast Surgery', 'Hematology and Oncology')
;


select department_specialty, note_type, count(distinct ir_id) as patient_count
from #TESTPATHOLNOTES
group by department_specialty, note_type
order by patient_count desc;

select *
from #TESTPATHOLNOTES;

-- 1246 patients in cohort 
