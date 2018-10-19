from django.conf.urls import url
from rest_framework.urlpatterns import format_suffix_patterns
from API import views

urlpatterns = [
    url(r'api/initFindings', views.initFindings),
    url(r'api/findings', views.findings),
    url(r'api/page', views.page),
    url(r'api/download', views.download)
    # url(r'api/connectDB', views.connectDB)
]