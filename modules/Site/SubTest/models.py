from django.db import models

# Create your models here.
class SubTestResult(models.Model):
    SCORES = [(i,i) for i in range(11)]
    TYPES = ((0,'128'),(1,'64'),(2,'d2'))
    DIST_TYPE = (('MS','Micro Shift'),('MR','Micro Rotate'),('BR','Blur'),('BL','Blocking'),('LCL','Light Change Linear'),('LCC','Light Change Const'))
    orgName = models.CharField(max_length=50)
    candName = models.CharField(max_length=50)
    distortion = models.CharField(max_length=2,choices = DIST_TYPE)
    #score = models.FloatField()
    score = models.IntegerField(choices=SCORES)
    imagetype = models.IntegerField(choices=TYPES)
    def __unicode__(self):
        return self.orgName + "_to_" + self.candName +"_with_"+self.distortion + "_"+ str(self.score)
