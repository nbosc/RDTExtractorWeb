#from snippets.models import Snippet
#from snippets.serializers import SnippetSerializer
from django.http import Http404
from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
import cx_Oracle
import pandas as pd
import json
import copy
from django.views.decorators.csrf import csrf_exempt
from .serializers import YourSerializer,FindingSerializer
from API.utils import extract

onto_df = pd.read_pickle("API/static/data/ontology.pkl")
compound_df = pd.read_pickle("API/static/data/compound.pkl")
findings_df = pd.read_pickle("API/static/data/findings.pkl.gz", compression='gzip')
study_df = pd.read_pickle("API/static/data/study.pkl")
merged_df = pd.merge(study_df[['study_id', 'subst_id', 'normalised_administration_route', 
                        'normalised_species', 'normalised_strain', 'exposure_period_days', 
                        'report_number']],
                    findings_df[['study_id', 'observation_normalised', 'grade',
                        'organ_normalised', 'dose', 'relevance', 'normalised_sex']],
                    how='left', on='study_id', left_index=False, right_index=False, 
                    sort=False)

@api_view(['GET'])
def source(request):

        global find_df

        host = '172.20.16.76'
        port = '1521'
        sid = 'ORA11G'
        user = 'vitic2016'
        password = 'T0Vitic2016'

        '''conn = connectDB(host, port, sid, user, password)

        cursor = conn.cursor()'''

        # Generate normalised study dataframe
        '''cmd = "SELECT DISTINCT PARENT_LUID AS study_id, RELEVANCE, \
            STANDARDISED_PATHOLOGY, STANDARDISED_ORGAN, DOSE \
            FROM HISTOPATHOLOGICALFI \
            WHERE STANDARDISED_PATHOLOGY IS NOT NULL \
            AND STANDARDISED_ORGAN IS NOT NULL"
        cursor.execute(cmd)
        results = cursor.fetchall()
        tmp_table = []
        for (study_id, relevance, observation, organ, dose) in results:
            tmp_table.append([study_id, relevance, observation, organ, dose])
            print (str(study_id)+","+str(relevance)+","+str(observation)+","+str(organ)+","+str(dose))
        find_df = pd.DataFrame(tmp_table, columns=['study_id', 'relevance',
                                                   'observation_normalised', 'organ_normalised', 'dose'])'''

        yourdata = [{"likes": 10, "comments": 0}, {"likes": 4, "comments": 23}]
        results = YourSerializer(yourdata, many=True).data
        return Response(yourdata)

@api_view(['GET'])
def initFindings(request):

    global merged_df

    page = int(request.GET.get("page"))

    num_studies = len(merged_df.study_id.unique().tolist())
    num_structures = len(merged_df.subst_id.unique().tolist())

    fullDict = {}
    fullDict['organs'] = merged_df.organ_normalised.dropna().unique().tolist()
    fullDict['organs'].sort()

    fullDict['observations'] = merged_df.observation_normalised.dropna().unique().tolist()
    fullDict['observations'].sort()

    fullDict['grade'] = merged_df.grade.dropna().unique().tolist()
    fullDict['grade'].sort()

    fullDict['routes'] = merged_df.normalised_administration_route.dropna().unique().tolist()
    fullDict['routes'].sort()

    fullDict['sex'] = merged_df.normalised_sex.dropna().unique().tolist()
    fullDict['sex'].sort()

    fullDict['species'] = merged_df.normalised_species.dropna().unique().tolist()
    fullDict['species'].sort()

    exposure_range = merged_df.exposure_period_days.dropna().unique().tolist()
    exposure_range.sort()
    fullDict['exposure_min'] = int(exposure_range[0])
    fullDict['exposure_max'] = int(exposure_range[-1])

    # Pagination
    if page != 0:
        init = (int(page) - 1) * 10
        end = init + 10
    else:
        init = 0
        end = len(merged_df)

    num_pages = int(len(merged_df) / 10)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = 7
    elif page > (num_pages - 4):
        previous = num_pages - 7
        nexts = num_pages
    else:
        previous = page - 4
        nexts = page + 3
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page >= num_pages):
        next_page = 0

    output = pd.merge(merged_df[init:end], compound_df[['subst_id','smiles']], how='left',
                        on='subst_id', left_index=False, right_index=False, sort=False)

    results = {
        'data': output.to_dict('records'),
        'allOptions': fullDict,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def findings(request):

    global merged_df, compound_df

    filtered_tmp = merged_df[:]

    # Sex
    sex = request.GET.getlist("sex")
    if sex:
        filtered_tmp = filtered_tmp[filtered_tmp.normalised_sex.str.lower() == sex[0].lower()]

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered_tmp = filtered_tmp[filtered_tmp.relevance == 'Treatment related']

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered_tmp = filtered_tmp[(filtered_tmp.exposure_period_days >= int(min_exposure)) &
                    (filtered_tmp.exposure_period_days <= int(max_exposure))]

    queryDict = {}
    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        queryDict['routes'] = 'normalised_administration_route == @all_routes'

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        queryDict['species'] = 'normalised_species == @all_species'

    # Organs
    all_organs = request.GET.getlist("organs")
    if len(all_organs) > 0:
        queryDict['organs'] = 'organ_normalised == @all_organs'
        
    # Observations
    all_observations = request.GET.getlist("observations")
    if len(all_observations) > 0:
        queryDict['observations'] = 'observation_normalised == @all_observations'

    # Grade
    all_grades = request.GET.getlist("grade")
    if len(all_grades) > 0:
        queryDict['grade'] = 'grade == @all_grades'

    #####################
    # Apply all filters #
    #####################
    query_string = ''
    if queryDict != {}:
        query_string = ' and '.join(list(queryDict.values()))
        filtered = filtered_tmp.query(query_string)
    else:
        filtered = filtered_tmp[:]

    num_studies = len(filtered.study_id.unique().tolist())
    num_structures = len(filtered.subst_id.unique().tolist())

    optionsDict = {}
    queryN = len(queryDict)
    if not filtered.empty:
        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('organs', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['organs'] = tmp_df.organ_normalised.dropna().unique().tolist()
        optionsDict['organs'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('observations', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['observations'] = tmp_df.observation_normalised.dropna().unique().tolist()
        optionsDict['observations'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('grade', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['grade'] = tmp_df.grade.dropna().unique().tolist()
        optionsDict['grade'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('routes', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['routes'] = tmp_df.normalised_administration_route.dropna().unique().tolist()
        optionsDict['routes'].sort()

        tmp_dict = copy.deepcopy(queryDict)
        tmp_dict.pop('species', None)
        valuesL = list(tmp_dict.values())
        if len(valuesL) > 0:
            query_string = ' and '.join(valuesL)
            tmp_df = filtered_tmp.query(query_string)
        else:
            tmp_df = filtered_tmp
        optionsDict['species'] = tmp_df.normalised_species.dropna().unique().tolist()
        optionsDict['species'].sort()

        exposure_range = filtered.exposure_period_days.dropna().unique().tolist()
        exposure_range.sort()
        optionsDict['exposure_min'] = int(exposure_range[0])
        optionsDict['exposure_max'] = int(exposure_range[-1])

        optionsDict['sex'] = merged_df.normalised_sex.dropna().unique().tolist()
        optionsDict['sex'].sort()

    ##############
    # Pagination #
    ##############

    page = int(request.GET.get("page"))

    if page != 0:
        init = (int(page) - 1) * 10
        end = init + 10
    else:
        init = 0
        end = len(filtered)

    num_pages = int(len(filtered) / 10)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = 7
    elif page > (num_pages - 4):

        previous = num_pages - 7
        nexts = num_pages
    else:
        previous = page - 4
        nexts = page + 3
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page >= num_pages):
        next_page = 0

    output = pd.merge(filtered[init:end], compound_df[['subst_id','smiles']], how='left', on='subst_id', left_index=False,
               right_index=False, sort=False)

    results = {
        'data': output.to_dict('records'),
        'allOptions': optionsDict,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def qualitative(request):


    global findings_df,study_df

    all_organs = request.GET.getlist("organs")
    all_observation = request.GET.getlist("observations")
    all_species = request.GET.getlist("species")
    all_routes = request.GET.getlist("routes")
    sex = request.GET.getlist("sex")
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    relevant = request.GET.get("treatmentRelated")

    #############
    # Filter    #
    #############

    filter_dict = {}
    if min_exposure:
        filter_dict['min_exposure'] = int(min_exposure[0])
    if max_exposure:
        filter_dict['max_exposure'] = int(max_exposure[0])
    if len(all_routes) > 0:
        filter_dict['route'] = all_routes
    if len(all_species) > 0:
        filter_dict['species'] = all_species


    filtered = extract.filter_study(filter_dict, study_df)
    filtered = findings_df[findings_df.study_id.isin(filtered.study_id)]

    filtered = pd.merge(filtered[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                  'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                               'normalised_species', 'exposure_period_days', 'report_number']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    if sex:
        filtered = filtered[filtered.normalised_sex.str.lower() == sex[0].lower()]
    if relevant:
        filtered = filtered[filtered.relevance == 'treatment related']
    if len(all_organs) > 0:
        filter_dict['organs'] = all_organs
        filtered = filtered[filtered['organ_normalised'].isin(all_organs)]
    if len(all_observation) > 0:
        filter_dict['observations'] = all_observation
        filtered = filtered[filtered['observation_normalised'].isin(all_observation)]


    ###################################
    # Get stats for relevant findings #
    ###################################

    # Get the number of studies per substance
    count_df = filtered.groupby(('subst_id')).study_id.nunique().to_frame().reset_index()
    # Get the global dose range per substance
    range_df = filtered[filtered.dose > 0]
    range_df = range_df.groupby(('subst_id')).dose.apply(extract.get_stats).unstack().reset_index()
    # Get all stats into a single dataframe
    stats_df = pd.merge(count_df, range_df, how='inner', on='subst_id',
                        left_index=False, right_index=False, sort=False)
    stats_df.columns = ['subst_id', 'study_count', 'dose_max', 'dose_min']

    # Expand
    filtered_piv = extract.expand(filtered, filter_dict, onto_df)

    ######################################
    # Aggragate by substance and finding #
    ######################################
    # Define finding as organ+observation
    filtered_piv.organ_normalised = filtered_piv.organ_normalised.fillna('NA')
    filtered_piv.observation_normalised = filtered_piv.observation_normalised.fillna('NA')
    filtered_piv['finding'] = filtered_piv.apply(
        lambda row: row.organ_normalised + '_' + row.observation_normalised, axis=1)
    filtered_piv = filtered_piv[['subst_id', 'finding', 'dose']]

    # Aggregate by substance and finding (as defined above), keeping the minimum dose
    # for each substance/finding instance
    group_df = filtered_piv.groupby(('subst_id', 'finding')).min().add_prefix('min_').reset_index()

    #######################################
    # Pivot so that each finding is a row #
    #######################################
    ### Qualitative
    group_df.loc[:, 'min_dose'] = 1
    pivotted_df = group_df.pivot_table(index='subst_id', columns='finding', values='min_dose').reset_index()
    pivotted_df['is_active'] = 'True'
    qualitative_df = pd.merge(stats_df, pivotted_df, how='left', on='subst_id',
                              left_index=False, right_index=False, sort=False)
    qualitative_df.is_active = qualitative_df.is_active.fillna('False')
    # Reorder columns
    cols = qualitative_df.columns.tolist()
    cols = cols[0:4] + [cols[-1]] + cols[4:-1]
    qualitative_df = qualitative_df[cols]

    print (qualitative_df.head())

    results = {
        "data" : qualitative_df.fillna(value=0)
    }


    return Response(results)

@api_view(['GET'])
def quantitative(request):


    global findings_df,study_df

    all_organs = request.GET.getlist("organs")
    all_observation = request.GET.getlist("observations")
    all_species = request.GET.getlist("species")
    all_routes = request.GET.getlist("routes")
    sex = request.GET.getlist("sex")
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    relevant = request.GET.get("treatmentRelated")

    global findings_df,study_df

    #############
    # Filter    #
    #############

    filter_dict = {}
    if min_exposure:
        filter_dict['min_exposure'] = int(min_exposure[0])
    if max_exposure:
        filter_dict['max_exposure'] = int(max_exposure[0])
    if len(all_routes) > 0:
        filter_dict['route'] = all_routes
    if len(all_species) > 0:
        filter_dict['species'] = all_species


    filtered = extract.filter_study(filter_dict, study_df)
    filtered = findings_df[findings_df.study_id.isin(filtered.study_id)]

    filtered = pd.merge(filtered[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                  'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                               'normalised_species', 'exposure_period_days', 'report_number']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    if sex:
        filtered = filtered[filtered.normalised_sex.str.lower() == sex[0].lower()]
    if relevant:
        filtered = filtered[filtered.relevance == 'treatment related']
    if len(all_organs) > 0:
        filter_dict['organs'] = all_organs
        filtered = filtered[filtered['organ_normalised'].isin(all_organs)]
    if len(all_observation) > 0:
        filter_dict['observations'] = all_observation
        filtered = filtered[filtered['observation_normalised'].isin(all_observation)]


    ###################################
    # Get stats for relevant findings #
    ###################################

    # Get the number of studies per substance
    count_df = filtered.groupby(('subst_id')).study_id.nunique().to_frame().reset_index()
    # Get the global dose range per substance
    range_df = filtered[filtered.dose > 0]
    range_df = range_df.groupby(('subst_id')).dose.apply(extract.get_stats).unstack().reset_index()
    # Get all stats into a single dataframe
    stats_df = pd.merge(count_df, range_df, how='inner', on='subst_id',
                        left_index=False, right_index=False, sort=False)
    stats_df.columns = ['subst_id', 'study_count', 'dose_max', 'dose_min']

    # Expand
    filtered_piv = extract.expand(filtered, filter_dict, onto_df)

    ######################################
    # Aggragate by substance and finding #
    ######################################
    # Define finding as organ+observation
    filtered_piv.organ_normalised = filtered_piv.organ_normalised.fillna('NA')
    filtered_piv.observation_normalised = filtered_piv.observation_normalised.fillna('NA')
    filtered_piv['finding'] = filtered_piv.apply(
        lambda row: row.organ_normalised + '_' + row.observation_normalised, axis=1)
    filtered_piv = filtered_piv[['subst_id', 'finding', 'dose']]

    # Aggregate by substance and finding (as defined above), keeping the minimum dose
    # for each substance/finding instance
    group_df = filtered_piv.groupby(('subst_id', 'finding')).min().add_prefix('min_').reset_index()

    #######################################
    # Pivot so that each finding is a row #
    #######################################
    ### Quantitative
    pivotted_df = group_df.pivot_table(index='subst_id', columns='finding', values='min_dose').reset_index()
    pivotted_df['is_active'] = 'True'
    quantitative_df = pd.merge(stats_df, pivotted_df, how='left', on='subst_id',
                               left_index=False, right_index=False, sort=False)
    quantitative_df.is_active = quantitative_df.is_active.fillna('False')
    # Reorder columns
    cols = quantitative_df.columns.tolist()
    cols = cols[0:4] + [cols[-1]] + cols[4:-1]
    quantitative_df = quantitative_df[cols]

    results = {
        "data": quantitative_df.fillna(value=0)
    }

    return Response(results)

@api_view(['GET'])
def study(request):

    id= request.GET.get("id")

    res = pd.merge(study_df[study_df.study_id == id], compound_df,
                    how='left', on='subst_id', left_index=False,
                    right_index=False, sort=False)

    results = {
        'study': res.fillna(value="-")
    }
    return Response(results)

# @api_view(['GET'])
# def connectDB(request):

#     host= request.GET.get("host")
#     port= int(request.GET.get("port"))
#     sid= request.GET.get("sid")
#     user= request.GET.get("user")
#     password= request.GET.get("password")

#     # try:
#     #     dsn_tns = cx_Oracle.makedsn(host, port, sid)
#     #     conn = cx_Oracle.connect(user, password, dsn=dsn_tns)        
#     # except:
#     #     #  cx_Oracle.DatabaseError as e:
#     #     conn = None
#     #     # error, = e.args
#     #     # connStatus = error.message
#     #     connStatus = 'Failed'
#     # else:
#     #     connStatus = 'Connected'

#     connStatus = 'Connected'

#     results = {
#         'connStatus': connStatus
#     }

#     send_data = FindingSerializer(results, many=False).data

#     return Response(send_data)