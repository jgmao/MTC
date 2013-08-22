from django import forms

class MatchForm(forms.Form):
    x =forms.IntegerField()
    y =forms.IntegerField()
    size = forms.IntegerField()