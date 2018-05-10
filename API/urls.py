from django.conf.urls import url
from rest_framework.urlpatterns import format_suffix_patterns
from API import views

urlpatterns = [
    url(r'findings', views.findings),
    url(r'source', views.source),
]