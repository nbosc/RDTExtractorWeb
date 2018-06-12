from django.conf.urls import url
from rest_framework.urlpatterns import format_suffix_patterns
from API import views

urlpatterns = [
    url(r'api/initFindings', views.initFindings),
    url(r'api/findings', views.findings),
    url(r'api/qualitative', views.qualitative),
    url(r'api/quantitative', views.quantitative),
    url(r'api/study', views.study)
    # url(r'api/connectDB', views.connectDB)
]