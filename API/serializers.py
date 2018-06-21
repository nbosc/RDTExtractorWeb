from rest_framework import serializers

class FindingSerializer (serializers.Serializer):
   data = serializers.ListField()
   allOptions = serializers.DictField()
   range_pages = serializers.ListField()
   num_pages = serializers.IntegerField()
   page = serializers.IntegerField()
   previous_page = serializers.IntegerField()
   next_page = serializers.IntegerField()
   num_studies = serializers.IntegerField()
   num_structures = serializers.IntegerField()
#    connStatus = serializers.CharField()