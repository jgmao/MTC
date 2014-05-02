#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dajaxice.decorators import dajaxice_register
from dajax.core import Dajax
from SubCoding.models import SubCodingResult
#from infant.forms import *
from django.core import serializers
@dajaxice_register
def say_hello(request,pk):
    dajax = Dajax()
    dajax.assign('#idUser', 'innerHTML', nap.user)
    dajax.assign('#idInfant', 'innerHTML', nap.infant)
    dajax.alert("Hello World!")
    return dajax.json()

@dajaxice_register
def multiply(request, a, b):
    dajax = Dajax()
    result = int(a) * int(b)
    dajax.assign('#result','value',str(result))
    return dajax.json()
