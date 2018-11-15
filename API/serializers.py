from rest_framework import serializers

class FindingSerializer (serializers.Serializer):

    data = serializers.ListField()
    allOptions = serializers.DictField()
    plotInfo = serializers.DictField()
    range_pages = serializers.ListField()
    num_pages = serializers.IntegerField()
    page = serializers.IntegerField()
    previous_page = serializers.IntegerField()
    next_page = serializers.IntegerField()
    num_studies = serializers.IntegerField()
    num_structures = serializers.IntegerField()
    num_findings = serializers.IntegerField()
    num_studies_positives =  serializers.IntegerField()
    num_structures_positives =  serializers.IntegerField()
    num_findings_positives =  serializers.IntegerField()
    num_studies_negatives =  serializers.IntegerField()
    num_structures_negatives =  serializers.IntegerField()
    num_findings_negatives =  serializers.IntegerField()
#    connStatus = serializers.CharField()

class InitFindingSerializer (serializers.Serializer):

    data = serializers.ListField()
    allOptions = serializers.DictField()
    plotInfo = serializers.DictField()
    range_pages = serializers.ListField()
    num_pages = serializers.IntegerField()
    page = serializers.IntegerField()
    previous_page = serializers.IntegerField()
    next_page = serializers.IntegerField()
    num_studies = serializers.IntegerField()
    num_structures = serializers.IntegerField()
    num_findings = serializers.IntegerField()

class Pageserializer (serializers.Serializer):

    data = serializers.ListField()
    range_pages = serializers.ListField()
    num_pages = serializers.IntegerField()
    page = serializers.IntegerField()
    previous_page = serializers.IntegerField()
    next_page = serializers.IntegerField()
#    connStatus = serializers.CharField()

