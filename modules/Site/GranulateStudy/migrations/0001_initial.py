# -*- coding: utf-8 -*-
import datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'Gran'
        db.create_table(u'GranulateStudy_gran', (
            (u'id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('orgName', self.gf('django.db.models.fields.CharField')(max_length=50)),
            ('gran', self.gf('django.db.models.fields.IntegerField')()),
        ))
        db.send_create_signal(u'GranulateStudy', ['Gran'])


    def backwards(self, orm):
        # Deleting model 'Gran'
        db.delete_table(u'GranulateStudy_gran')


    models = {
        u'GranulateStudy.gran': {
            'Meta': {'object_name': 'Gran'},
            'gran': ('django.db.models.fields.IntegerField', [], {}),
            u'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'orgName': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        }
    }

    complete_apps = ['GranulateStudy']