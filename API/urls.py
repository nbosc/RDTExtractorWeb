from django.conf.urls import url
from rest_framework.urlpatterns import format_suffix_patterns
from API import views

urlpatterns = [
    url(r'api/findings/organs', views.organs),
    url(r'api/findings/observations', views.observations),
    url(r'api/findings/species', views.species),
    url(r'api/findings/routes', views.routes),
    url(r'api/findings/sex', views.sex),
    url(r'api/findings', views.findings),

]