from django import forms
from django.forms import ModelForm
from models import SubTestResult
from django.forms.fields import DateField, ChoiceField, MultipleChoiceField
from django.forms.widgets import RadioSelect, CheckboxSelectMultiple
#class SubTestForm(ModelForm):
#    class Meta:
#        model =  SubTestResult
#        fields = ['score']

class SubTestForm(forms.Form):
    score = forms.ChoiceField(widget=forms.RadioSelect, choices=SubTestResult.SCORES)