# parity_parsenotes
import re
# import numpy as np
# import pandas as pd

def clean_string(notestring):
    """
    This function takes the raw string of the genetic test result section as input
    and cleans it by:
    - removing any remnants of RTF formats (this is no longer necessary after I refined my SQL query)
    - changing all sequences of whitespaces to a single space 
    - stripping whitespace at the edges
    """
    # cleaned_string = re.sub('\\\\par+', '', notestring)
    # cleaned_string = re.sub('\\\\plain\\\\f\d\\\\fs22|\\\\lang1033\\\\hich\\\\f\d\\\\dbch\\\\f\d\\\\loch\\\\f\d\\\\cf\d\\\\fs22|\\\\[ib]\d*', '', cleaned_string)
    cleaned_string = re.sub('\s+', ' ', notestring) # replace all whitespaces with single space
    cleaned_string = cleaned_string.strip() # strip outer whitespace 
    return cleaned_string

def strip_lower(string):
    return re.sub('\s', '', string).lower()


def get_gravidapara(notestring, resultdict):
    gp0 = re.findall(r'gravida ?[0-9]+[,\? ]*para ?[0-9]', notestring, re.IGNORECASE) # ,?\?? ?
    gp1 = re.findall(r'g:? ?[0-9],?/? ?p:? ?[0-9]', notestring, re.IGNORECASE)
    gp2 = re.findall(r'p ?[0-9],?/? ?g ?[0-9]', notestring, re.IGNORECASE)
    gp3 = re.findall(r'Grav[ida]* Para Term.* [0-9] [0-9]', notestring) # Pt is a OB History 
    gp4 = re.findall(r'Grav[ida]* Para Term Preterm Ab[ortions]* TAB SAB Ect[opic]* Mult[iple]* Living [0-9] [a-z]', notestring, re.IGNORECASE)
    gp5 = re.findall(r'nulliparous|nuliparious|nullparous|nuliparous|Pt has not been pregnant|Pt has never been pregnant', notestring, re.IGNORECASE)
    gp6 = re.findall(r'Grav[ida]* Para Term Preterm Ab[ortions]* TAB SAB Ect[opic]* Mult[iple]* Living [a-z]', notestring, re.IGNORECASE)
    gp7 = re.findall(r'Pt is a OB History \?|Pt is a OB History No data available|Pt is a OB History No obstetric history on file|Pt is a OB History None|has not been asked about pregnancy', notestring, re.IGNORECASE)
    if len(gp0)>0: # gravida, para
        gplist = re.findall(r'[0-9]+', gp0[0])
        resultdict['gravida'] = int(gplist[0])
        resultdict['para'] = int(gplist[1])
    elif len(gp1)>0: # G, P
        gplist = re.findall(r'[0-9]+', gp1[0])
        resultdict['gravida'] = int(gplist[0])
        resultdict['para'] = int(gplist[1])
    elif len(gp2)>0: # P, G
        gplist = re.findall(r'[0-9]+', gp2[0])
        resultdict['gravida'] = int(gplist[1])
        resultdict['para'] = int(gplist[0])
    elif len(gp3)>0: # Gravida Para Term... <num> <num>
        gplist = re.findall(r'[0-9]+', gp3[0])
        resultdict['gravida'] = int(gplist[0])
        resultdict['para'] = int(gplist[1])
    elif len(gp4)>0: # Gravida Para Term... <num>
        gplist = re.findall(r'[0-9]+', gp4[0])
        resultdict['gravida'] = int(gplist[0])
        resultdict['para'] = int(gplist[0])
    elif len(gp5)>0: # nulliparous 
        resultdict['gravida'] = 0
        resultdict['para'] = 0
    elif (len(gp6)>0 | len(gp7)>0): # Gravida Para Term... Pt is... -or- OB History missing 
        resultdict['gravida'] = None
        resultdict['para'] = None
    # else:
    #     resultdict['gravida'] = None
    #     resultdict['para'] = None
    return resultdict


def get_menarche(notestring, resultdict):
    menarche = re.findall(r'menarchea?t?a?g?e?:?.?[0-9]+|ageofmenarche[0-9]+|ageatmenarche[0-9]+', strip_lower(notestring))
    if len(menarche)!=0: # if menarche not found
        resultdict['menarche'] = int(re.findall(r'[0-9]+', menarche[0])[0])
    else:
        menarche = re.findall(r'menarchea?t?a?g?e?:?\?|menarchea?t?a?g?e?:?\*+|menarcheatanunknownage|menarcheinherteens', strip_lower(notestring))
        if len(menarche)>0:
            resultdict['menarche'] = 'unknown'
        else:
            menarche = re.findall(r'menarcheyoungerthan[0-9]+|menarchebefore[0-9]+|menarchepriortoage[0-9]+', strip_lower(notestring))
            if len(menarche)>0:
                resultdict['menarche'] = int(re.findall(r'[0-9]+', menarche[0])[0])
    return resultdict


def get_menopause(notestring, resultdict):
    premenopausal = re.findall(r'ageofmenopause:?n/?a|ageofmenopause:?premenopausal|pre\-?menopausal|menstruatesregularly|regularmenstration|stillmenstu?rating|regularmensturation', strip_lower(notestring))
    postmenopausal = re.findall(r'ageofmenopause:?[0-9]+|menopauseatage[0-9]+|menopauseat[0-9]+|menopauseatageof[0-9]+|menopausein[0-9]+|lastperiodwasatage[0-9]+|ageatmenopause:[0-9]+|ispost\-?menopausal', strip_lower(notestring))
    # peri = re.findall(r'peri\-?menopausal', strip_lower(notestring))
    # hyst = re.findall(r'hysterectomy|hystrectomy', notestring, re.IGNORECASE)
    if len(premenopausal)>0:
        resultdict['menopause'] = 'premenopausal'
    # elif len(peri)>0:
    #     resultdict['menopause'] = 'perimenopausal'
    # elif len(hyst)>0:
    #     resultdict['menopause'] = 'hysterectomy'
    elif len(postmenopausal)>0:
        postmenopausal = re.findall(r'[0-9]+', postmenopausal[0])
        if len(postmenopausal)>0:
            if int(postmenopausal[0])<30:
                resultdict['menopause'] = None
            elif int(postmenopausal[0])>1000:
                resultdict['menopause'] = 'in' + postmenopausal[0]
            elif int(postmenopausal[0])>=100:
                resultdict['menopause'] = int(postmenopausal[0][:2])
            else:
                resultdict['menopause'] = int(postmenopausal[0])
        else:
            resultdict['menopause'] = 'postmenopausal'
    else:
        # Consultation 2
        premenopausal = re.findall(r'Menstruates: Regularly: \[x\]', notestring, re.IGNORECASE)
        postmenopausal = re.findall(r'Post-menopausal \[x\]', notestring, re.IGNORECASE)
        perimenopausal = re.findall(r'Perimenopausal: \[x\]', notestring, re.IGNORECASE)
        if len(premenopausal)>0:
            resultdict['menopause'] = 'premenopausal'
        if len(postmenopausal)>0:
            resultdict['menopause'] = 'postmenopausal'
        if len(perimenopausal)>0:
            resultdict['menopause'] = 'perimenopausal'
    return resultdict


def get_lmp(notestring, resultdict):
    lmp = re.findall(r'lastmenstrualperiodwas[0-9\/\-]+|lmp:[0-9\/\-]+|lastperiodwason[0-9\/\-]+|lmpwas[0-9\/\-]+|lastperiodwas[0-9\/\-]+|lmpon[0-9\/\-]+', strip_lower(notestring))
    if len(lmp)!=0:
        resultdict['lmp'] = re.findall(r'[0-9\/\-]+', lmp[0])[0]
    else:
        lmp = re.findall(r'No LMP recorded', notestring, re.IGNORECASE)
        if len(lmp)>0:
            resultdict['lmp'] = 'notrecorded'
        # else:
        #     resultdict['lmp'] = None
    return resultdict


def get_ageoffirstpreg(notestring, resultdict):
    agefirstpreg = re.findall(r'[0-9]*[na\/]*yearsoldattimeoffirstpregnancy|firstfulltermpregancyattheageof[0-9]+|firstfulltermpregnancyattheageof[0-9]+|firstfulltermpregnancyattheageof[na\/]|firstlivechildbirthatage[0-9]+|ageoffirstlivebirth:?[0-9]+', strip_lower(notestring))
    if len(agefirstpreg)>0:
        agefirstpreg = re.findall(r'[0-9]+',agefirstpreg[0])
        if len(agefirstpreg)>0:
            resultdict['agefirstpreg'] = int(agefirstpreg[0])
        # else:
        #     resultdict['agefirstpreg'] = None
    else:
        agefirstpreg = re.findall(r'First live birth - age [0-9]+', notestring, re.IGNORECASE)
        if len(agefirstpreg)>0:
            resultdict['agefirstpreg'] = int(re.findall(r'[0-9]+', agefirstpreg[0])[0])
        # else:
        #     resultdict['agefirstpreg'] = None
    return resultdict


def get_ageoflastpreg(notestring, resultdict):
    agelastpreg = re.findall(r'Age of last live birth: [0-9]+', notestring, re.IGNORECASE)
    if len(agelastpreg)>0:
        agelastpreg = int(re.findall(r'[0-9]+', agelastpreg[0])[0])
        resultdict['ageoflastpreg'] = agelastpreg
    return resultdict


def get_lacthist(notestring, resultdict):
    # for Initial Consultation
    lacthist = re.findall(r'Lactation history:', notestring) # .*breastcancer
    if len(lacthist)!=0:
        lacthist = re.split(r'Lactation history:|Breast Cancer:', notestring)[1].strip()
        if re.match(r'^none|^n\/a|^no|^denies|^never', lacthist, re.IGNORECASE):
            resultdict['lastnursed'] = 'none'
        else:
            ## First, find pattern like 'xx years ago'
            lastnursed = re.findall(r'[0-9\.\/ ]+ years ago|[0-9\.\/ ]+ years prior|time since last nursed [0-9\.\/ ]+ years', lacthist)
            if len(lastnursed)>0:
                try:
                    resultdict['lastnursed'] = '%syearsago' % int(re.findall(r'[0-9\.]', lacthist)[0])
                except:
                    if re.findall(r'[0-9\.\/ ]', lacthist)[0]=='1 1/2':
                        resultdict['lastnursed'] = '1.5yearsago'
                    else:
                        print('Unrecognized pattern: %s' % lacthist)
            else:
                ## Next, find patterns like 'xx months ago'
                # print(lacthist, '\n')
                lastnursed = re.findall(r'[0-9]+ months ago', lacthist)
                if len(lastnursed)!=0:
                    resultdict['lastnursed'] = '%smonthsago' % int(re.findall(r'[0-9]+', lastnursed[0])[0])
                else:
                    ## Finally, find patterns where the year is shown 
                    lastnursed = re.findall(r'in [0-9][0-9][0-9][0-9]|as [0-9][0-9][0-9][0-9]', notestring)
                    if len(lastnursed)!=0:
                        resultdict['lastnursed'] = 'in%s' % int(re.findall(r'[0-9][0-9][0-9][0-9]', lastnursed[0])[0])
                    else:
                        resultdict['lastnursed'] = 'missing'
    else:
        # for Consultation 1
        lacthist = re.findall(r'cumulative nursing history of [0-9]+ months? and last nursed [0-9\<\.~]+ years? ago', notestring, re.IGNORECASE)
        if len(lacthist)==0:
            lacthist = re.findall(r'cumulative nursing history of [0-9]+ months? and last nursed in [0-9]+', notestring, re.IGNORECASE)
            if len(lacthist)==0:
                lacthist = re.findall(r'cumulative nursing history of [0-9]+ months?', notestring, re.IGNORECASE)
                if len(lacthist)!=0:
                    resultdict['duration'] = int(re.findall(r'[0-9]+', lacthist[0])[0])
            else:
                resultdict['duration'] = '%smonths' % int(re.findall(r'[0-9]+', lacthist[0])[0])
                resultdict['lastnursed'] = 'in' + re.findall(r'[0-9]+', lacthist[0])[1] # concatenate so that it looks like 'in2005' or something like this
        else:
            resultdict['duration'] = '%smonths' % int(re.findall(r'[0-9]+', lacthist[0])[0])
            resultdict['lastnursed'] = '%syearsago' % re.findall(r'[0-9]+', lacthist[0])[1]
    return resultdict


def read_initcons(ir_id, title_type, notestring):
    resultdict = {}
    resultdict['ir_id'] = ir_id
    resultdict['title_type'] = title_type
    ## Get gravida/para
    resultdict = get_gravidapara(notestring, resultdict)
    if 'gravida' not in resultdict.keys():
        if re.match(r'ptis[as]*obhistorygravparatermpretermabortions', strip_lower(notestring)): # 1.3
            # the "Pt is a s OB" appears in only one case and is nulliparous
            resultdict['gravida'] = 0
            resultdict['para'] = 0
        elif re.match(r'ptisa*g[0-9]', strip_lower(notestring)): # 1.6 - only gravida available
            g = re.findall(r'ptisa*g[0-9]', strip_lower(notestring))[0]
            resultdict['gravida'] = int(g[-1])
        elif re.match(r'ptisagravida[0-9],parida[0-9]', strip_lower(notestring)): # 1.9
            # This is a weird case
            resultdict['gravida'] = 2
            resultdict['para'] = 2
    ## Get menarche 
    resultdict = get_menarche(notestring, resultdict)
    ## Get menopause 
    resultdict = get_menopause(notestring, resultdict)
    ## get age of first preg 
    resultdict = get_ageoffirstpreg(notestring, resultdict)
    ## last nursed 
    resultdict = get_lacthist(notestring, resultdict)
    ## LMP
    resultdict = get_lmp(notestring, resultdict)
    # last birth
    # looks like no last birth info available in this format
    return resultdict

def read_con_newpat(ir_id, title_type, notestring):
    resultdict = {}
    resultdict['ir_id'] = ir_id
    resultdict['title_type'] = title_type
    ## Get gravida/para
    resultdict = get_gravidapara(notestring, resultdict)
    ## Get menarche 
    resultdict = get_menarche(notestring, resultdict)
    ## Get menopause 
    resultdict = get_menopause(notestring, resultdict)
    ## get age of first preg 
    resultdict = get_ageoffirstpreg(notestring, resultdict)
    ## last nursed 
    resultdict = get_lacthist(notestring, resultdict)
    ## LMP
    resultdict = get_lmp(notestring, resultdict)
    # last birth
    resultdict = get_ageoflastpreg(notestring, resultdict)
    return resultdict
