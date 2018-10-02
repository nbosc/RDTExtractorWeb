# from snippets.models import Snippet
# from snippets.serializers import SnippetSerializer
from django.http import Http404
from rest_framework.response import Response
from rest_framework.decorators import api_view
from rest_framework import status
import pandas as pd
# Disable SettingWithCopyWarning warnings
pd.set_option('chained_assignment', None)
import json
import copy
import math
import pickle
import time, cProfile
from django.views.decorators.csrf import csrf_exempt
from .serializers import FindingSerializer,Pageserializer
from API.utils import extract

def get_stats(group):
    return {'min': group.min(), 'max': group.max()}

# Load dataframes with information on studies, compounds and findings
t0 = time.time()
substance_df = pd.read_pickle("API/static/data/substance.pkl")
# t1 = time.time()
study_df = pd.read_pickle("API/static/data/study.pkl")
# t2 = time.time()
organ_onto_df = pd.read_pickle("API/static/data/organ_ontology.pkl")
# t3 = time.time()
observation_onto_df = pd.read_pickle("API/static/data/observation_ontology.pkl")
# t4 = time.time()
all_df = pd.read_pickle("API/static/data/all.pkl.gz", compression="gzip")
tf = time.time()
print ('Loading:\n\t{}'.format(tf-t0))
# print ('substance only:\n\t{}'.format(t1-t0))
# print ('study only:\n\t{}'.format(t2-t1))
# print ('organ only:\n\t{}'.format(t3-t2))
# print ('observation only:\n\t{}'.format(t4-t3))
# print ('all only:\n\t{}'.format(tf-t4))

@api_view(['GET'])
def initFindings(request):

    global output_df, optionsDict

    t0 = time.time()
    output_df = pd.read_pickle("API/static/data/output.pkl")
    results = pickle.load(open("API/static/data/init_results.pkl", 'rb'))
    optionsDict = results['allOptions']
    send_data = FindingSerializer(results, many=False).data
    tf = time.time()
    print ('init:\n\t{}'.format(tf-t0))

    return Response(send_data)

@api_view(['GET'])
def findings(request):

    global all_df, substance_df, study_df, output_df, optionsDict
    t0 = time.time()

    #####################
    # Apply all filters #
    #####################

    ##
    ## Substance-level filters
    ##
    filtered_subs = substance_df[:]

    # Pharmacological action
    all_pharm = request.GET.getlist("pharmacological_action")
    if len(all_pharm) > 0:
        filtered_subs = filtered_subs[filtered_subs.targetAction.isin(all_pharm)]
    t1 = time.time()

    # Compound name
    all_compound_name = request.GET.getlist("compound_name")
    if len(all_compound_name) > 0:
        # Solve issue with plus signs in compound names being converted to spaces
        plus_signs = [x.replace(' ', '+') for x in all_compound_name]
        all_compound_name = all_compound_name+plus_signs
        all_compound_name = list(set(all_compound_name))
        filtered_subs = filtered_subs[filtered_subs.common_name.isin(all_compound_name)]

    # CAS number
    all_cas_number = request.GET.getlist("cas_number")
    if len(all_cas_number) > 0:
        filtered_subs = filtered_subs[filtered_subs.cas_number.isin(all_cas_number)]

    ##
    ## Study-level filters
    ##
    filtered_studies = study_df[:]
    filtered_studies = filtered_studies[filtered_studies.subst_id.isin(filtered_subs.subst_id)]
    t2 = time.time()

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered_studies = filtered_studies[(filtered_studies.exposure_period_days >= int(min_exposure)) &
                    (filtered_studies.exposure_period_days <= int(max_exposure))]
    t3 = time.time()

    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        filtered_studies = filtered_studies[filtered_studies.admin_route.isin(all_routes)]
    t4 = time.time()

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        filtered_studies = filtered_studies[filtered_studies.species.isin(all_species)]
    t5 = time.time()

    ##
    ## Finding-level filters
    ##
    filtered = pd.merge(all_df, filtered_studies[['study_id']], on='study_id', how='inner')
    t6 = time.time()

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered = filtered[filtered.relevance == 'Treatment related']
    t7 = time.time()

    # Sex
    sex = request.GET.get("sex")
    if sex:
        filtered = filtered[filtered.sex == sex]
    t8 = time.time()

    ##
    ## Filter parameters, observations and grades by category
    ##
    # Use this df to store each parameter / observation and add them together
    # so they don't become mutually exclusive
    additive_df = pd.DataFrame(columns=filtered.columns)
    
    # parameters and observations
    all_parameters = request.GET.getlist("parameters")
    tmp_parameters_dict = {}
    all_categories = set([])
    if len(all_parameters) > 0:
        all_parameters = all_parameters[0].split('@')
        for v in all_parameters:
            category, val = v.split('|')
            # Expand based on the ontology
            expanded_val = list(organ_onto_df[organ_onto_df.parent_term == val].child_term)
            expanded_val = [val]+expanded_val
            if category not in tmp_parameters_dict:
                tmp_parameters_dict[category] = expanded_val
            else:
                tmp_parameters_dict[category].extend(expanded_val)
        all_categories = set(tmp_parameters_dict.keys())
    t9 = time.time()

    # Observations
    all_observations = request.GET.getlist("observations")
    tmp_observations_dict = {}
    if len(all_observations) > 0:
        all_observations = all_observations[0].split('@')
        for v in all_observations:
            category, val = v.split('|')
            category = category.strip()
            val = val.strip()
            # Expand based on the ontology
            expanded_val = list(observation_onto_df[observation_onto_df.parent_term == val].child_term)
            expanded_val = [val]+expanded_val
            if category not in tmp_observations_dict:
                tmp_observations_dict[category] = expanded_val
            else:
                tmp_observations_dict[category].extend(expanded_val)
            if category not in tmp_observations_dict:
                tmp_observations_dict[category] = [val]
            else:
                tmp_observations_dict[category].append(val)
        all_categories = all_categories.union(set(tmp_parameters_dict.keys()))
    t10 = time.time()

    if bool(all_categories):
        # At least one observation or parameter filter has been applied
        for category in all_categories:
            if category in tmp_parameters_dict and category in tmp_observations_dict:
                or_df = filtered[(filtered.endpoint_type == category.strip()) &
                                (filtered.parameter.isin(tmp_parameters_dict[category])) &
                                (filtered.observation.isin(tmp_observations_dict[category]))]
            elif category in tmp_parameters_dict:
                or_df = filtered[(filtered.endpoint_type == category.strip()) &
                                (filtered.parameter.isin(tmp_parameters_dict[category]))]
            elif category in tmp_observations_dict:
                or_df = filtered[(filtered.endpoint_type == category) & 
                                (filtered.observation.isin(tmp_observations_dict[category]))]
            additive_df = pd.concat([additive_df, or_df])
        filtered = additive_df[:]
    t11 = time.time()

    #############
    # Aggregate #
    #############

    num_studies = filtered.study_id.nunique()
    num_structures = filtered.subst_id.nunique()
    t12a = time.time()

    filtered_findings = filtered[['dose', 'observation', 'parameter', 'relevance', 'sex', 'endpoint_type', 'study_id']].drop_duplicates().groupby(['dose', 'observation', 'parameter', 'relevance', 'sex', 'endpoint_type', 'study_id'])
    t12b = time.time()
    num_findings = filtered_findings.ngroups
    t12c = time.time()
    filtered_findings = filtered_findings.count().reset_index()
    t12d = time.time()
    t12 = time.time()

    #############
    # Plot info #
    #############

    plot_info = {}
    # Species
    species = filtered.groupby(['species','study_id']).count()
    species = species.reset_index().groupby(['species'])['species'].count()
    species.sort_values(ascending=False, inplace=True)
    plot_info['species'] = [[],[]]
    sum_value = 0
    for index, value in species.iteritems():
        p = float(value)/num_studies
        if p < 0.015:
            plot_info['species'][0].append('Other')
            plot_info['species'][1].append(num_studies-sum_value)
            break
        sum_value += value
        plot_info['species'][0].append(index)
        plot_info['species'][1].append(value)

    # Treatment related
    relevance = filtered_findings.groupby(['relevance'])['relevance'].count()
    relevance.sort_values(ascending=False, inplace=True)
    plot_info['relevance'] = [relevance.index, relevance.values]

    # Source
    source = filtered_findings.groupby(['endpoint_type'])['endpoint_type'].count()
    source.sort_values(ascending=False, inplace=True)
    plot_info['source'] = [[], []]
    sum_value = 0
    for index, value in source.iteritems():
        p = float(value)/num_findings
        if p < 0.015:
            plot_info['source'][0].append('Other')
            plot_info['source'][1].append(num_findings-sum_value)
            break
        sum_value += value
        plot_info['source'][0].append(index)
        plot_info['source'][1].append(value)
    t13 = time.time()

    #################
    # Create output #
    #################

    if not filtered.empty:
        study_count_df = filtered.dropna(subset=['species'])[['subst_id', 'species', 'study_id']].groupby(['subst_id', 'species']).study_id.nunique().reset_index()
        study_count_df.columns = ['subst_id', 'species', 'study_count']
        study_count_df.loc[:,'count'] = study_count_df.species + ': ' + study_count_df.study_count.astype(int).astype(str)
        study_count_df = study_count_df[['subst_id', 'count']].groupby('subst_id').agg(lambda x : '\n'.join(x)).reset_index()

        output_df = pd.merge(filtered[['subst_id']].drop_duplicates(),
                        substance_df[['subst_id', 'cas_number', 'common_name',  
                                'smiles', 'status', 'targetActionList']].drop_duplicates(), 
                        how='left', on='subst_id', left_index=False, right_index=False, 
                        sort=False)
        output_df = pd.merge(output_df, study_count_df,
                        how='left', on='subst_id', left_index=False, right_index=False, 
                        sort=False)
        output_df = output_df.drop_duplicates()
        output_df.common_name = output_df.common_name.str.replace(', ', '\n')
        t14 = time.time()
    else:
        output_df = pd.DataFrame(columns=['subst_id', 'cas_number', 'common_name', 
                                    'smiles', 'status', 'targetActionList', 'count'])

    ##############
    # Pagination #
    ##############

    page = int(request.GET.get("page"))

    if page != 0:
        init = (page - 1) * 5
        end = init + 5
    else:
        init = 0
        end = len(filtered)

    num_pages = math.ceil(output_df.subst_id.nunique() / 5.)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = min([7, num_pages])
    elif page > (num_pages - 4):
        previous = num_pages - 7
        nexts = num_pages
    else:
        previous = page - 4
        nexts = page + 3
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page > num_pages):
        next_page = 0

    results = {
        'data': output_df[init:end].fillna(value="-").to_dict('records'),
        'allOptions': optionsDict,
        'plotInfo': plot_info,
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page,
        'num_studies': num_studies,
        'num_structures': num_structures,
        'num_findings' : num_findings
    }

    send_data = FindingSerializer(results, many=False).data
    tf = time.time()
    print ('TOTAL: %.4f' %(tf-t0))    
    print ('\tpharm action: %.4f' %(t1-t0))
    print ('\tstudies filter: %.4f' %(t2-t1))
    print ('\texposure: %.4f' %(t3-t2))
    print ('\troutes: %.4f' %(t4-t3))
    print ('\tspecies: %.4f' %(t5-t4))
    print ('\tfindings filer: %.4f' %(t6-t5))
    print ('\trelevant: %.4f' %(t7-t6))
    print ('\tsex: %.4f' %(t8-t7))
    print ('\tparamenters: %.4f' %(t9-t8))
    print ('\tobservations: %.4f' %(t10-t9))
    print ('\tapply parameters and observations: %.4f' %(t11-t10))
    print ('\taggregate: %.4f' %(t12-t11))
    print ('\t\tnunique: %.4f' %(t12a-t11))
    print ('\t\tgroupby: %.4f' %(t12b-t12a))
    print ('\t\tngroups: %.4f' %(t12c-t12b))
    print ('\t\tcount reset index: %.4f' %(t12d-t12c))
    print ('\tplot info: %.4f' %(t13-t12))
    print ('\tfinal step: %.4f' %(t14-t13))
    return Response(send_data)

@api_view(['GET'])
def page(request):

    global output_df

    page = int(request.GET.get("page"))

    total = output_df.subst_id.nunique()
    if page != 0:
        init = (page - 1) * 5
        end = init + 5
    else:
        init = 0
        end = total

    num_pages = math.ceil(total / 5.)

    # Range of pages to show
    if page < 4:
        previous = 0
        nexts = min([7, num_pages])
    elif page > (num_pages - 4):
        previous = num_pages - 7
        nexts = num_pages
    else:
        previous = page - 4
        nexts = page + 3
    range_pages = range(1, num_pages + 1)[previous:nexts]
    previous_page = page - 1
    next_page = page + 1

    if (next_page > num_pages):
        next_page = 0

    results = {
        'data': output_df[init:end].fillna(value="-").to_dict('records'),
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': page,
        'previous_page': previous_page,
        'next_page': next_page
    }

    send_data = Pageserializer(results, many=False).data
    return Response(send_data)
