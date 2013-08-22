# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):

        # Changing field 'SubTestResult.score'
        db.alter_column(u'SubTest_subtestresult', 'score', self.gf('django.db.models.fields.IntegerField')())

    def backwards(self, orm):

        # Changing field 'SubTestResult.score'
        db.alter_column(u'SubTest_subtestresult', 'score', self.gf('django.db.models.fields.FloatField')())

    models = {
        u'SubTest.subtestresult': {
            'Meta': {'object_name': 'SubTestResult'},
            'candName': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'distortion': ('django.db.models.fields.CharField', [], {'max_length': '2'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'orgName': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'score': ('django.db.models.fields.IntegerField', [], {})
        }
    }

    complete_apps = ['SubTest']