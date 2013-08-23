#!/usr/bin/env python
# -*- coding: utf-8 -*-

from django import forms
from django.forms import ModelForm
from models import Gran
from django.forms.fields import DateField, ChoiceField, MultipleChoiceField
from django.forms.widgets import RadioSelect, CheckboxSelectMultiple

class GranTestForm(forms.Form):
    gran = forms.ChoiceField(widget=forms.RadioSelect, choices=Gran.SIZE)
