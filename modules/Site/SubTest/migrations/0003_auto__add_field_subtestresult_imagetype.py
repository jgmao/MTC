# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding field 'SubTestResult.imagetype'
        db.add_column(u'SubTest_subtestresult', 'imagetype',
                      self.gf('django.db.models.fields.IntegerField')(default=0),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting field 'SubTestResult.imagetype'
        db.delete_column(u'SubTest_subtestresult', 'imagetype')


    models = {
        u'SubTest.subtestresult': {
            'Meta': {'object_name': 'SubTestResult'},
            'candName': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'distortion': ('django.db.models.fields.CharField', [], {'max_length': '2'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'imagetype': ('django.db.models.fields.IntegerField', [], {}),
            'orgName': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'score': ('django.db.models.fields.IntegerField', [], {})
        }
    }

    complete_apps = ['SubTest']