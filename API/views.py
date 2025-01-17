import os
import math
import pickle
import io
import zipfile
import time, cProfile
import pandas as pd
# Disable SettingWithCopyWarning warnings
pd.set_option('chained_assignment', None)
from .serializers import FindingSerializer, Pageserializer, InitFindingSerializer
from django.http import HttpResponse
from rest_framework.response import Response
from rest_framework.decorators import api_view

def get_stats(group):
    return {'min': group.min(), 'max': group.max()}

caption = ''

# Load dataframes with information on studies, compounds and findings
substance_df = pd.read_pickle("API/static/data/substance.pkl")
study_df = pd.read_pickle("API/static/data/study.pkl")
study_cmpd_df = pd.read_pickle("API/static/data/study_cmpd.pkl")
organ_onto_df = pd.read_pickle("API/static/data/organ_ontology.pkl")
observation_onto_df = pd.read_pickle("API/static/data/observation_ontology.pkl")
all_df = pd.read_pickle("API/static/data/all.pkl.gz", compression="gzip")
filtered = all_df[:]

optionsDict = {}

@api_view(['GET'])
def initFindings(request):

    global output_df, filtered, optionsDict

    output_df = pd.read_pickle("API/static/data/output.pkl")
    filtered = all_df[:]
    results = pickle.load(open("API/static/data/init_results.pkl", 'rb'))
    optionsDict = results['allOptions']
    send_data = FindingSerializer(results, many=False).data

    return Response(send_data)

@api_view(['GET'])
def findings(request):

    global all_df, substance_df, study_df, study_cmpd_df, output_df, optionsDict, filtered, caption

    #####################
    # Apply all filters #
    #####################

    ##
    ## Substance-level filters
    ##
    filtered_subs = substance_df[:]
    substance_caption = ''

    # Pharmacological action
    all_pharm = request.GET.getlist("pharmacological_action")
    if len(all_pharm) > 0:
        filtered_subs = filtered_subs[filtered_subs.targetAction.isin(all_pharm)]
        substance_caption += '\tPharmacological action: %s\n' %', '.join(all_pharm)

    # Compound name
    all_compound_name = request.GET.getlist("compound_name")
    if len(all_compound_name) > 0:
        # Solve issue with plus signs in compound names being converted to spaces
        plus_signs = [x.replace(' ', '+') for x in all_compound_name]
        all_compound_name = all_compound_name+plus_signs
        all_compound_name = list(set(all_compound_name))
        filtered_subs = filtered_subs[filtered_subs.common_name.isin(all_compound_name)]
        substance_caption += '\tCompound name: %s\n' %', '.join(all_compound_name)

    # CAS number
    all_cas_number = request.GET.getlist("cas_number")
    if len(all_cas_number) > 0:
        filtered_subs = filtered_subs[filtered_subs.cas_number.isin(all_cas_number)]
        substance_caption += '\tCas number: %s\n' %', '.join(all_cas_number)

    ##
    ## Study-level filters
    ##
    filtered_studies = study_df[:]
    # Create filtered negative studies
    filtered_studies_negative = study_cmpd_df[:]

    filtered_studies = filtered_studies[filtered_studies.nonStandard_inchi_key.isin(filtered_subs.nonStandard_inchi_key)]
    study_caption = ''

    # Species
    all_species = request.GET.getlist("species")
    if len(all_species) > 0:
        filtered_studies = filtered_studies[filtered_studies.species.isin(all_species)]
        filtered_studies_negative = filtered_studies_negative[filtered_studies_negative.species.isin(all_species)]
        study_caption += '\tSpecies: %s\n' %', '.join(all_species)

    # Administration route
    all_routes = request.GET.getlist("routes")
    if len(all_routes) > 0:
        filtered_studies = filtered_studies[filtered_studies.admin_route.isin(all_routes)]
        filtered_studies_negative = filtered_studies_negative[filtered_studies_negative.admin_route.isin(all_routes)]
        study_caption += '\tAdminitration route: %s\n' %', '.join(all_routes)

    # Exposure
    min_exposure = request.GET.get("min_exposure")
    max_exposure = request.GET.get("max_exposure")
    if min_exposure and max_exposure:
        filtered_studies = filtered_studies[(filtered_studies.exposure_period_days >= int(min_exposure)) &
                    (filtered_studies.exposure_period_days <= int(max_exposure))]
        study_caption += '\tExposure range: %d to %d days\n' %(int(min_exposure), int(max_exposure))

    # Dose
    min_dose = request.GET.get("min_dose")
    max_dose = request.GET.get("max_dose")
    if min_dose and max_dose:
        filtered_studies = filtered_studies[(filtered_studies.dose_min >= float(min_dose)) &
                    (filtered_studies.dose_max <= int(max_dose))]
        study_caption += '\tDose range: %s to %s mg/kg\n' %(min_dose, max_dose)        
    # Min negative Dose
    min_negative_dose = request.GET.get("negative_min_dose")
    if min_negative_dose:
        study_caption += '\tMinimum tested dose for negatives: %s mg/kg\n' %(min_negative_dose)

    ##
    ## Finding-level filters
    ##
    filtered = pd.merge(all_df, filtered_studies[['study_id']], on='study_id', how='inner')
    if min_negative_dose:
        filtered_negative = pd.merge(all_df, filtered_studies_negative[['study_id']], 
                                    on='study_id', how='inner')
    finding_caption = ''

    # Sex
    sex = request.GET.get("sex")
    if sex:
        filtered = filtered[filtered.sex == sex]
        if sex == 'F':
            finding_caption += '\tSex: Female\n'
        elif sex == 'M':
            finding_caption += '\tSex: Male\n'
        else:
            finding_caption += '\tSex: Both\n'

    # Relevancy
    relevant = request.GET.get("treatmentRelated")
    if relevant:
        filtered = filtered[filtered.relevance]
        finding_caption += '\tTreatment-related only\n'

    #
    # Filter parameters, observations and grades by category
    #
    # Use this df to store each parameter / observation and add them together
    # so they don't become mutually exclusive
    filtered['positive'] = True
    additive_df = pd.DataFrame(columns=filtered.columns)
    additive_df.positive = additive_df.positive.astype('bool')
    finding_observation_caption = ''
    
    # Parameters
    all_parameters = request.GET.getlist("parameters")
    tmp_parameters_dict = {}
    all_categories = set([])
    if len(all_parameters) > 0:
        all_parameters = all_parameters[0].split('@')
        finding_observation_caption += '\tParameter: %s\n' %', '.join(all_parameters)
        for v in all_parameters:
            category, val = v.split('|')
            category = category.strip()
            all_categories.add(category)
            val = val.strip()
            checked = []
            expanded = set([])
            to_check = [val]
            # Expand recursively based on the ontology
            while to_check:
                val = to_check.pop()
                expanded.add(val)
                checked.append(val)
                children = set(organ_onto_df[organ_onto_df.parent_term == val].child_term)
                expanded.update(children)
                for child in children:
                    if child not in checked:
                        to_check.append(child)
            if category not in tmp_parameters_dict:
                tmp_parameters_dict[category] = list(expanded)
            else:
                tmp_parameters_dict[category].extend(list(expanded))

    # Observations
    all_observations = request.GET.getlist("observations")
    tmp_observations_dict = {}
    if len(all_observations) > 0:
        all_observations = all_observations[0].split('@')
        finding_observation_caption += '\tObservation: %s\n' %', '.join(all_observations)
        for v in all_observations:
            category, val = v.split('|')
            category = category.strip()
            all_categories.add(category)
            val = val.strip()
            checked = []
            expanded = set([])
            to_check = [val]
            # Expand recursively based on the ontology
            while to_check:
                val = to_check.pop()
                expanded.add(val)
                checked.append(val)
                children = set(observation_onto_df[observation_onto_df.parent_term == val].child_term)
                expanded.update(children)
                for child in children:
                    if child not in checked:
                        to_check.append(child)
            if category not in tmp_observations_dict:
                tmp_observations_dict[category] = list(expanded)
            else:
                tmp_observations_dict[category].extend(list(expanded))

    if bool(all_categories):
        # At least one observation or parameter filter has been applied
        for category in all_categories:
            category = category.strip()
            if category in tmp_parameters_dict and category in tmp_observations_dict:
                or_df = filtered[(filtered.endpoint_type == category) &
                                (filtered.parameter.isin(tmp_parameters_dict[category])) &
                                (filtered.observation.isin(tmp_observations_dict[category])) &
                                (filtered.dose > 0)]
            elif category in tmp_parameters_dict:
                or_df = filtered[(filtered.endpoint_type == category) &
                                (filtered.parameter.isin(tmp_parameters_dict[category])) &
                                (filtered.dose > 0)]
            elif category in tmp_observations_dict:
                or_df = filtered[(filtered.endpoint_type == category) &
                                (filtered.observation.isin(tmp_observations_dict[category])) &
                                (filtered.dose > 0)]

            or_df['positive'] = True
            additive_df = pd.concat([additive_df, or_df])

            if min_negative_dose:
                tmp_findings_df = all_df[(all_df.endpoint_type == category) &
                                       (all_df.parameter.isin(tmp_parameters_dict[category])) &
                                       (all_df.dose > 0)]
                tmp_study_cmpd_df = filtered_studies_negative[(filtered_studies_negative.study_id.isin(tmp_findings_df.study_id))]
                negative_df = filtered_studies_negative[~(filtered_studies_negative.nonStandard_inchi_key.isin(tmp_study_cmpd_df.nonStandard_inchi_key)) &(filtered_studies_negative.dose_max >= float(min_negative_dose))]
                negative_df = all_df[all_df.nonStandard_inchi_key.isin(negative_df.nonStandard_inchi_key)]
                negative_df['positive'] = False
                additive_df = pd.concat([additive_df, negative_df])

        filtered = additive_df[:]

    ##################################
    # Generate caption summarizing   #
    # the filtering criteria applied #
    ##################################

    caption = ''
    if study_caption != '':
        caption += 'Study-level filters:\n%s\n' %study_caption
    if finding_observation_caption != '' or finding_caption != '':
        caption += 'Finding-level filters:\n'
        if finding_observation_caption != '':
            caption += '%s' %finding_observation_caption
        if finding_caption != '':
            caption += '%s' %finding_caption
        caption += '\n'
    if substance_caption != '':
        caption += 'Substance-level filters:\n%s\n' %substance_caption
    
    #############
    # Aggregate #
    #############

    num_studies = filtered.study_id.nunique()
    num_structures = filtered.nonStandard_inchi_key.nunique()

    num_studies_positives = filtered.study_id[filtered.positive].nunique()
    num_structures_positives = filtered.nonStandard_inchi_key[filtered.positive].nunique()
    num_studies_negatives = filtered.study_id[~filtered.positive].nunique()
    num_structures_negatives = filtered.nonStandard_inchi_key[~filtered.positive].nunique()

    filtered_findings = filtered[['dose', 'observation', 'parameter', 'relevance', 'sex',
                            'endpoint_type', 'study_id', 'positive']].drop_duplicates().groupby(['dose', 'observation', 'parameter', 'relevance', 'sex', 'endpoint_type', 'study_id', 'positive'])
    num_findings = filtered_findings.ngroups

    filtered_findings = filtered_findings.count().reset_index()
    num_findings_positive= len(filtered_findings[filtered_findings.positive])
    num_findings_negatives = len(filtered_findings[~filtered_findings.positive])

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
    plot_info['relevance'] = [[], []]
    for index, value in relevance.iteritems():
        if index:
            plot_info['relevance'][0].append('Treatment-related')
        else:
            plot_info['relevance'][0].append('Not related')
        plot_info['relevance'][1].append(value)

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

    #################
    # Create output #
    #################

    if not filtered.empty:

        study_count_df = filtered.dropna(subset=['species'])[['nonStandard_inchi_key', 'species', 'study_id']].groupby(['nonStandard_inchi_key', 'species']).study_id.nunique().reset_index()
        study_count_df.columns = ['nonStandard_inchi_key', 'species', 'study_count']
        study_count_df.loc[:,'count'] = study_count_df.species + ': ' + study_count_df.study_count.astype(int).astype(str)
        study_count_df = study_count_df[['nonStandard_inchi_key', 'count']].groupby('nonStandard_inchi_key').agg(lambda x : '\n'.join(x)).reset_index()

        output_df = pd.merge(filtered[['nonStandard_inchi_key','positive']].drop_duplicates(),
                        substance_df[['nonStandard_inchi_key', 'cas_number', 'common_name',  
                                'smiles', 'status', 'targetActionList']].drop_duplicates(), 
                        how='left', on='nonStandard_inchi_key', left_index=False, right_index=False, 
                        sort=False)
        output_df = pd.merge(output_df, study_count_df,
                        how='left', on='nonStandard_inchi_key', left_index=False, right_index=False, 
                        sort=False)
        output_df = output_df.drop_duplicates()
        output_df.common_name = output_df.common_name.str.replace(', ', '\n')
    else:
        output_df = pd.DataFrame(columns=['nonStandard_inchi_key', 'cas_number', 'common_name', 
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

    num_pages = math.ceil(output_df.nonStandard_inchi_key.nunique() / 5.)

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

    output_df.to_pickle("../output_df.pkl")
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
        'num_findings' : num_findings,
        'num_studies_positives': num_studies_positives,
        'num_structures_positives': num_structures_positives,
        'num_findings_positives': num_findings_positive,
        'num_studies_negatives': num_studies_negatives,
        'num_structures_negatives': num_structures_negatives,
        'num_findings_negatives': num_findings_negatives
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)

@api_view(['GET'])
def page(request):

    global output_df

    page = int(request.GET.get("page"))

    total = output_df.nonStandard_inchi_key.nunique()
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

@api_view(['GET'])
def download(request):

    global substance_df, filtered, caption

    ####################################
    ## Generate output files          ##
    ## the filtering criteria applied ##
    ####################################

    smiles_df = substance_df[:]
    smiles_df = smiles_df[['nonStandard_inchi_key', 'std_smiles']].drop_duplicates()
    ids_df = substance_df[['nonStandard_inchi_key', 'subst_id']].drop_duplicates()
    ids_df = ids_df.groupby(['nonStandard_inchi_key'],as_index=False).agg(lambda x: ', '.join(x))
    positives = filtered.groupby('nonStandard_inchi_key')['positive'].first().to_frame().reset_index()
    output_df = filtered[:]

    # Define finding as organ+observation
    output_df[output_df.positive].dropna(subset=['parameter', 'observation'], inplace=True)
    output_df['finding'] = output_df.parameter+'_'+output_df.observation
    quant_filtered_df = output_df[['nonStandard_inchi_key', 'finding', 'dose','positive']]
    quant_filtered_df = quant_filtered_df[quant_filtered_df.dose>0]
    
    ##
    ## Get stats for relevant findings
    ##
    group_df = output_df.groupby(('nonStandard_inchi_key'))
    # Get the number of studies per substance
    count_df = group_df.study_id.nunique().to_frame().reset_index()
    min_df = group_df.dose_min.min().to_frame().reset_index()
    max_df = group_df.dose_max.max().to_frame().reset_index()
    group_df = quant_filtered_df[quant_filtered_df.positive].groupby(('nonStandard_inchi_key'))
    min_observation_dose_df = group_df.dose.min().to_frame().reset_index()


    # Get all stats into a single dataframe
    stats_df = pd.merge(count_df, min_df, how='inner', on='nonStandard_inchi_key', 
                        left_index=False, right_index=False, sort=False)
    stats_df = pd.merge(stats_df, max_df, how='inner', on='nonStandard_inchi_key', 
                        left_index=False, right_index=False, sort=False)
    stats_df = pd.merge(stats_df, min_observation_dose_df, 
                        how='left', on='nonStandard_inchi_key', 
                        left_index=False, right_index=False, sort=False)
    stats_df.columns = ['nonStandard_inchi_key', 'study_count', 'dose_min', 
                        'dose_max', 'min_observation_dose']
    
    ##
    ## Aggragate by substance and finding
    ##

    # Aggregate by substance and finding (as defined above), 
    # keeping the minimum dose for each substance/finding instance
    group_df = quant_filtered_df[quant_filtered_df.dose > 0].groupby(('nonStandard_inchi_key', 'finding')).min().add_prefix('min_').reset_index()
    
    ##
    ## Pivot so that each finding is a row
    ##
    
    ### Quantitative
    pivotted_df = group_df[group_df.min_positive].pivot_table(index='nonStandard_inchi_key', columns='finding',
                                                           values='min_dose').reset_index()
    negative = pd.DataFrame(quant_filtered_df[~quant_filtered_df.positive]['nonStandard_inchi_key'].unique(), columns=['nonStandard_inchi_key'])
    pivotted_df = pd.concat([pivotted_df, negative], ignore_index=True)

    quantitative_df = pd.merge(stats_df, pivotted_df, how='left', on='nonStandard_inchi_key', 
                                left_index=False, right_index=False, sort=False)
    # Reorder columns
    cols = quantitative_df.columns.tolist()
    quantitative_df = quantitative_df[cols]
    quantitative_df = pd.merge(quantitative_df, ids_df, 
                               how='left', on='nonStandard_inchi_key', 
                               left_index=False, right_index=False, sort=False)
    quantitative_df = pd.merge(quantitative_df, 
                               smiles_df[['nonStandard_inchi_key', 'std_smiles']],
                               how='left', on='nonStandard_inchi_key', 
                               left_index=False, right_index=False, sort=False)
    quantitative_df = pd.merge(quantitative_df, positives,
                               how='left', on='nonStandard_inchi_key',
                               left_index=False, right_index=False, sort=False)

    ### Qualitative
    group_df = output_df.groupby(['nonStandard_inchi_key', 'finding','positive']).study_id.nunique().reset_index(name='counts')
    pivotted_df = group_df[group_df.positive].pivot_table(index='nonStandard_inchi_key', columns='finding', values='counts').reset_index()
    pivotted_df = pd.concat([pivotted_df, negative], ignore_index=True)
    qualitative_df = pd.merge(stats_df, pivotted_df, how='left', on='nonStandard_inchi_key',
                                left_index=False, right_index=False, sort=False)
    # Reorder columns
    cols = qualitative_df.columns.tolist()
    cols = cols[0:5]+[cols[-1]]+cols[5:-1]
    qualitative_df = qualitative_df[cols]
    qualitative_df = pd.merge(qualitative_df, ids_df, how='left', on='nonStandard_inchi_key', 
                               left_index=False, right_index=False, sort=False)
    qualitative_df = pd.merge(qualitative_df, 
                              smiles_df[['nonStandard_inchi_key', 'std_smiles']],
                              how='left', on='nonStandard_inchi_key', 
                              left_index=False, right_index=False, sort=False)
    qualitative_df = pd.merge(qualitative_df, positives,
                              how='left', on='nonStandard_inchi_key',
                              left_index=False, right_index=False, sort=False)
    
    ##
    ## Create the HttpResponse object with the appropriate CSV header.
    ##
    files = {}

    buffer = io.StringIO()
    qualitative_df.to_csv(buffer, encoding='utf-8', sep='\t', index=False)
    buffer.seek(0)
    files['qualitative'] = buffer

    buffer = io.StringIO()
    quantitative_df.to_csv(buffer, encoding='utf-8', sep='\t', index=False)
    buffer.seek(0)
    files['quantitative'] = buffer

    zipped_file = io.BytesIO()
    with zipfile.ZipFile(zipped_file, "a", zipfile.ZIP_DEFLATED, False) as zipper:
        for i, file in files.items():
            file.seek(0)
            file_text = caption+file.read()
            zipper.writestr("{}.tsv".format(i), file_text)
    zipped_file.seek(0)

    response = HttpResponse(zipped_file, content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="results.zip"'

    return response