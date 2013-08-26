from django.db import models

# Create your models here.
class Gran(models.Model):
  orgName = models.CharField(max_length=50)
  SIZE = ((0,4),(1,8),(2,16),(3,32),(4,64),(5,128),(6,-1))#the last one means none of any size
  create_date = models.DateTimeField(auto_now_add=True)
  gran = models.IntegerField(choices=SIZE)
  ip = models.CharField(max_length=50)
  def __unicode__(self):
    return self.orgName+"_"+str(self.gran)

