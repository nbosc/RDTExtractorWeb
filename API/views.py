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

findings_df = pd.read_pickle("API/static/data/findings.pkl.gz", compression='gzip')
ontology_df = pd.read_pickle("API/static/data/ontology.pkl")
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

    # getting values from post
    page = request.GET.get('page')
    page=0
    # Just first time
    global filtered

    if page == 0:
        filtered = pd.merge(findings_df[['study_id', 'observation_normalised', 'organ_normalised', 'dose', \
                                      'relevance', 'normalised_sex']],
                            study_df[['study_id', 'subst_id', 'normalised_administration_route', \
                                   'normalised_species', 'exposure_period_days']],
                            how='left', on='study_id', left_index=False,
                            right_index=False, sort=False)
        page = 1

    #Range of page
    init = (int(page) - 1) * 10;
    end = init + 10

    # adding the values in a context variable
    num_pages = int(len(filtered) / 10)

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

        'data': filtered[0:10].to_dict('records'),
        'range_pages': range_pages,
        'num_pages': num_pages,
        'page': 1,
        'previous_page': previous_page,
        'next_page': next_page,
    }

    send_data = FindingSerializer(results, many=False).data
    return Response(send_data)


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