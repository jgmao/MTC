# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'SubTestResult'
        db.create_table(u'SubTest_subtestresult', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('orgName', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('candName', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('distortion', self.gf('django.db.models.fields.CharField')(max_length=2)),
            ('score', self.gf('django.db.models.fields.FloatField')()),
        ))
        db.send_create_signal(u'SubTest', ['SubTestResult'])


    def backwards(self, orm):
        # Deleting model 'SubTestResult'
        db.delete_table(u'SubTest_subtestresult')


    models = {
        u'SubTest.subtestresult': {
            'Meta': {'object_name': 'SubTestResult'},
            'candName': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'distortion': ('django.db.models.fields.CharField', [], {'max_length': '2'}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'orgName': ('django.db.models.fields.CharField', [], {'max_length': '50'}),
            'score': ('django.db.models.fields.FloatField', [], {})
        }
    }

    complete_apps = ['SubTest']