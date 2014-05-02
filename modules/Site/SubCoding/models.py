from django.db import models

# Create your models here.

class SubCodingResult(models.Model):
    SCORES = [(i,i) for i in range(11)]
    SIZES = ((0,'128'),(1,'64'),(2,'32'),(3,'16'))
    orgName = models.CharField(max_length=50)
    candName = models.CharField(max_length=50)
    score = models.IntegerField(choices=SCORES)
    create_date = models.DateTimeField(auto_now_add=True)
    ip = models.CharField(max_length=50)
    def __unicode__(self):
        return self.orgName + "_to_" + self.candName +"_"+ str(self.score)
