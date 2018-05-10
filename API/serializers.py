from rest_framework import serializers

class YourSerializer(serializers.Serializer):
   """Your data serializer, define your fields here."""
   comments = serializers.IntegerField()
   likes = serializers.IntegerField()

class FindingSerializer (serializers.Serializer):

   data = serializers.ListField(True)
   range_pages = serializers.ListField()
   num_pages = serializers.IntegerField()
   page = serializers.IntegerField()
   previous_page = serializers.IntegerField()
   next_page = serializers.IntegerField()