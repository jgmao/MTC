from django.conf.urls import patterns, include, url
from django.conf import settings
from django.conf.urls.static import static
# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'TPCSite.views.home', name='home'),
    # url(r'^TPCSite/', include('TPCSite.TPCSite.urls')),
    #url(r'^static\/(?P<path>.*)$','django.views.static.serve'),
    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
    #url(r'^$','Matching.views.show_image',name='show_image'),
    #url(r'^$','Matching.views.show_result',name='result'),
    url(r'^matching/','Matching.views.matching',{'x':None,'y':None,'size':None,},name='matching'),
    url(r'^result/','Matching.views.show_result',name='result'),
    #url(r'^start64/','SubTest.views.startpage64',name='teststart64'),
    #url(r'^startd2/','SubTest.views.startpaged2',name='teststartd2'),
    #url(r'^starti64/','SubTest.views.startinter',name='teststarti64'),
    #url(r'^start128/','SubTest.views.startpage',name='teststart128'),
    url(r'^subtest/(\d{1,2})/$','SubTest.views.show_image',name='subtest'),
    #url(r'^threadtest/$','SubTest.views.threadtest',name='threadtest'),
    url(r'^output/$','SubTest.views.output',name='output'),
    url(r'^granout/$','GranulateStudy.views.output', name='granout'),
    url(r'^dist/(\w+)/(\d{2,3})/(\w+)/$','SubTest.views.startcompare',name='startcompare'),
    url(r'^granstart/$','GranulateStudy.views.initTestData',name='startgran'),
    url(r'^grantest/$','GranulateStudy.views.show_image',name='grantest'),
    url(r'^subcoding/','SubCoding.views.startsession',name='subcoding'),
    url(r'^codingprocess/','SubCoding.views.codingprocess',name='codingprocess'),
    url(r'^cand_selected/','SubCoding.views.cand_selected',name='cand_selected'),
)
urlpatterns += staticfiles_urlpatterns()
urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
