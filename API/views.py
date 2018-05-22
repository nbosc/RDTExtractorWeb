#from snippets.models import Snippet
#from snippets.serializers import SnippetSerializer
from django.http import Http404
from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
import cx_Oracle
import pandas as pd
import json
from django.views.decorators.csrf import csrf_exempt
from .serializers import YourSerializer,FindingSerializer
from API.utils import extract

findings_df = pd.read_pickle("API/static/data/findings.pkl.gz", compression='gzip')
onto_df = pd.read_pickle("API/static/data/ontology.pkl")
study_df = pd.read_pickle("API/static/data/study.pkl")


@api_view(['GET'])
def source(request):

        global find_df

        host = '172.20.16.76'
        port = '1521'
        sid = 'ORA11G'
        user = 'vitic2016'
        password = 'T0Vitic2016'

        conn = connectDB(host, port, sid, user, password)

        cursor = conn.cursor()

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
def findings(request):


    global findings_df,study_df

    all_organs = request.GET.getlist("organs")
    all_observation = request.GET.getlist("observations")
    all_species = request.GET.getlist("species")
    all_routes = request.GET.getlist("routes")
    sex = request.GET.getlist("sex")
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    relevant = request.GET.get("treatmentRelated")
    page = int(request.GET.get("page"))

    global findings_df,study_df

    #Filter
    ########################################3
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

    num_studies = len(filtered.study_id.unique().tolist())
    num_structures = len(filtered.subst_id.unique().tolist())

    # Range of page

    if page != 0:
        init = (int(page) - 1) * 20;
        end = init + 20
    else:
        init = 0
        end = len(filtered)
    # adding the values in a context variable
    num_pages = int(len(filtered) / 20)

    # paginator = Paginator(df.to_dict('records'), 10)


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

    results = {

        'data': filtered[init:end].to_dict('records'),
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
def organs(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    organs=filtered.organ_normalised.unique().tolist()
    #organs.remove(None)
    organs.sort()

    results = {
        'organs': organs
    }
    return Response(results)

@api_view(['GET'])
def observations(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    observations = filtered.observation_normalised.unique().tolist()
    # observations.remove(None)
    observations.sort()

    results = {
        'observations': observations
    }
    return Response(results)

@api_view(['GET'])
def routes(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    route = filtered.normalised_administration_route.unique().tolist()
    route.remove(None)
    route.sort()

    results = {
        'route': route
    }

    return Response(results)

@api_view(['GET'])
def species(request):
    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    species = filtered.normalised_species.unique().tolist()
    species.remove(None)
    species.sort()

    results = {
        'species': species
    }

    return Response(results)\

@api_view(['GET'])
def sex(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    sex = filtered.normalised_sex.unique().tolist()
    sex.remove(None)
    sex.sort()

    results = {
        'sex': sex
    }

    return Response(results)

@api_view(['GET'])
def organs(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    organs=filtered.organ_normalised.unique().tolist()
    #organs.remove(None)
    organs.sort()

    results = {
        'organs': organs
    }
    return Response(results)


@api_view(['GET'])
def observations(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    observations = filtered.observation_normalised.unique().tolist()
    # observations.remove(None)
    observations.sort()

    results = {
        'observations': observations
    }
    return Response(results)

@api_view(['GET'])
def routes(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    route = filtered.normalised_administration_route.unique().tolist()
    route.remove(None)
    route.sort()

    results = {
        'route': route
    }

    return Response(results)

@api_view(['GET'])
def species(request):
    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    species = filtered.normalised_species.unique().tolist()
    species.remove(None)
    species.sort()

    results = {
        'species': species
    }

    return Response(results)\

@api_view(['GET'])
def sex(request):

    global findings_df, study_df

    filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                     'relevance', 'normalised_sex']],
                        study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                  'normalised_species', 'exposure_period_days']],
                        how='left', on='study_id', left_index=False,
                        right_index=False, sort=False)

    sex = filtered.normalised_sex.unique().tolist()
    sex.remove(None)
    sex.sort()

    results = {
        'sex': sex
    }

    return Response(results)

def testConnectDB (host, port, sid, user, password) :
    conn= None
    dsn_tns = cx_Oracle.makedsn(host,port,sid)
    try:
        conn = cx_Oracle.connect(user,password,dsn=dsn_tns)
    except cx_Oracle.DatabaseError as e:
        result=e
    finally:
        if conn is not None:
            conn.close()
            result="Database connected."
    return result

def connectDB (host, port, sid, user, password):

    conn = None
    dsn_tns = cx_Oracle.makedsn(host, port, sid)
    try:
        conn = cx_Oracle.connect(user, password, dsn=dsn_tns)
    except cx_Oracle.DatabaseError as e:
        return None
    return conn